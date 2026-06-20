"""
    ccr2d(decay_expts, buildup_expts, T)

Start interactive GUI for analyzing 2D measurements of cross-correlated relaxation data.

Fits intensities of peaks in a series of 2D spectra to a model of the form:

    tanh(η * T) = I_buildup / I_decay

or

    tanh(η * T) = sqrt((I_buildup1 * I_buildup2) / (I_decay1 * I_decay2))


# Arguments
- `decay_expts`: List of NMR data files for decay experiments
- `buildup_expts`: List of NMR data files for buildup experiments
- `T`: Relaxation time constant (in seconds)

# Example:
```julia
ccr2d("decay_expt", "buildup_expt", 0.08)    # single decay and buildup experiment

ccr2d(["decay_expt1", "decay_expt2"],        # symmetric reconversion experiments
      ["buildup_expt1", "buildup_expt2"], 0.08)
```
"""
function ccr2d(decay_expts::Vector, buildup_expts::Vector, T)
    expt = CCRExperiment(decay_expts, buildup_expts, T)
    return gui!(expt)
end

ccr2d(decay_expt::String, buildup_expt::String, T) = ccr2d([decay_expt], [buildup_expt], T)

"""
    CCRExperiment

Cross-correlated relaxation experiment with decay and buildup spectra.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `isbuildup`: Vector of booleans indicating buildup (true) or decay (false) spectra
- `T`: Relaxation time constant in seconds
- `issymmetric`: Whether symmetric reconversion is used (pairs of experiments)
"""
struct CCRExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    isbuildup::Any
    T::Float64
    issymmetric::Bool

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function CCRExperiment(specdata, peaks, isbuildup, T, issymmetric)
        expt = new(specdata, peaks, isbuildup, T, issymmetric,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}())
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

struct CCRVisualisation <: VisualisationStrategy end
visualisationtype(::CCRExperiment) = CCRVisualisation()
primaryparam(::CCRExperiment) = :eta

"""
    CCRExperiment(decay_expts, buildup_expts, T)

Create CCR experiment from lists of decay and buildup experiment files.

# Arguments
- `decay_expts`: Vector of file paths to decay experiments
- `buildup_expts`: Vector of file paths to buildup experiments
- `T`: Relaxation time constant in seconds

For symmetric reconversion, provide pairs of experiments (2 decay + 2 buildup).
For standard CCR, provide single experiments (1 decay + 1 buildup).
"""
function CCRExperiment(decay_expts::Vector, buildup_expts::Vector, T)
    # Validate input
    length(decay_expts) == length(buildup_expts) ||
        throw(ArgumentError("Number of decay and buildup experiments must match"))
    length(decay_expts) in (1, 2) ||
        throw(ArgumentError("CCR requires either 1 or 2 pairs of decay/buildup experiments"))

    issymmetric = length(decay_expts) == 2

    # Merge experiments: interleave decay and buildup
    planefilenames = String[]
    isbuildup = Bool[]
    for i in eachindex(decay_expts)
        push!(planefilenames, decay_expts[i])
        push!(isbuildup, false)
        push!(planefilenames, buildup_expts[i])
        push!(isbuildup, true)
    end

    specdata = preparespecdata(planefilenames, isbuildup, CCRExperiment)
    peaks = Observable(Vector{Peak}())

    return CCRExperiment(specdata, peaks, isbuildup, T, issymmetric)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(planefilenames, isbuildup, ::Type{CCRExperiment})
    @debug "Preparing spec data for CCR experiment: $planefilenames"
    spectra = map(loadnmr, planefilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(bu -> bu ? "Buildup" : "Decay", isbuildup)

    return SpecData(spectra, x, y,
                    z ./ σ[1],
                    σ ./ σ[1],
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::CCRExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.0)
    R2y0 = MaybeVector(10.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    newpeak.postparameters[:eta] = Parameter("CCR rate eta", 0.0)
    newpeak.postparameters[:amp] = Parameter("Reference amplitude", maximum(amp0))

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::CCRExperiment, xbounds=nothing, ybounds=nothing)
    n = length(z)
    for i in 1:n
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds[i]]
        y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds[i]]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]
        R2x = peak.parameters[:R2x].value[][i]
        R2y = peak.parameters[:R2y].value[][i]
        amp = peak.parameters[:amp].value[][i]
        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(2π * hz(x0, xaxis), R2x, 2π * hz(xs, xaxis),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x * R2y) *
             NMRTools.NMRBase._lineshape(2π * hz(y0, yaxis), R2y, 2π * hz(ys, yaxis),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting.

Computes CCR rate η from the model:
- Single experiments: tanh(η * T) = I_buildup / I_decay
- Symmetric (pairs): tanh(η * T) = sqrt((I_buildup1 * I_buildup2) / (I_decay1 * I_decay2))
"""
function postfit!(peak::Peak, expt::CCRExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog = 10

    # Get intensities with uncertainties
    A = peak.parameters[:amp].value[] .± peak.parameters[:amp].uncertainty[]

    # Separate buildup and decay intensities
    Idecay = A[expt.isbuildup .== false]
    Ibuildup = A[expt.isbuildup .== true]

    # Calculate ratio depending on symmetric or standard CCR
    if expt.issymmetric
        # Symmetric reconversion: sqrt((I_buildup1 * I_buildup2) / (I_decay1 * I_decay2))
        ratio = sqrt((Ibuildup[1] * Ibuildup[2]) / (Idecay[1] * Idecay[2]))
    else
        # Standard CCR: I_buildup / I_decay
        ratio = Ibuildup[1] / Idecay[1]
    end

    # Calculate eta: tanh(η * T) = ratio → η = atanh(ratio) / T
    eta = atanh(ratio) / expt.T

    # Store results
    peak.postparameters[:eta].value[] .= Measurements.value(eta)
    peak.postparameters[:eta].uncertainty[] .= Measurements.uncertainty(eta)

    # Reference amplitude is the mean decay intensity
    Iref = sum(Idecay) / length(Idecay)
    peak.postparameters[:amp].value[] .= Measurements.value(Iref)
    peak.postparameters[:amp].uncertainty[] .= Measurements.uncertainty(Iref)

    return peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::CCRExperiment, idx)
    if expt.isbuildup[idx]
        "Buildup ($idx of $(nslices(expt)))"
    else
        "Decay ($idx of $(nslices(expt)))"
    end
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::CCRExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "η: $(peak.postparameters[:eta].value[][1] ± peak.postparameters[:eta].uncertainty[][1]) s⁻¹\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "Amplitude: $(peak.postparameters[:amp].value[][1] ± peak.postparameters[:amp].uncertainty[][1])\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::CCRExperiment)
    mode = expt.issymmetric ? "Symmetric reconversion" : "Standard"
    return "Analysis type: CCR ($mode)\n" *
           "Relaxation time T: $(expt.T) s\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Experiments: $(join(expt.specdata.zlabels, ", "))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

## visualisation
function get_ccr_data(peak, expt::CCRExperiment)
    x = 1:nslices(expt)

    if isnothing(peak)
        y = fill(0.0, nslices(expt))
        err = [(1.0 * i, 0.0, 0.0) for i in 1:nslices(expt)]
        # Return default fit lines (won't be visible at y=0 when bars are at 0)
        fit = [0.0, 0.0]
        return (x, y, err, fit)
    end

    # Normalize to decay amplitude
    Iref = peak.postparameters[:amp].value[][1]

    y = [peak.parameters[:amp].value[][i] / Iref for i in 1:nslices(expt)]
    err = [(1.0 * i,
            peak.parameters[:amp].value[][i] / Iref,
            peak.parameters[:amp].uncertainty[][i] / Iref)
           for i in 1:nslices(expt)]

    # Calculate expected buildup/decay ratio from fitted eta
    eta = peak.postparameters[:eta].value[][1]
    expected_ratio = tanh(eta * expt.T)
    # fit line shows: decay = 1.0, buildup = expected_ratio
    fit = [1.0, expected_ratio]

    return (x, y, err, fit)
end

function completestate!(state, expt, ::CCRVisualisation)
    @debug "completing state for CCR visualisation"
    state[:peak_plot_data] = lift(peak -> get_ccr_data(peak, expt), state[:current_peak])
    state[:peak_plot_x] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_y] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_err] = lift(d -> d[3], state[:peak_plot_data])
    return state[:peak_plot_fit] = lift(d -> d[4], state[:peak_plot_data])
end

function plot_peak!(panel, peak, expt, ::CCRVisualisation)
    x, y, err, fit = get_ccr_data(peak, expt)

    # Color bars: blue for decay, orange for buildup
    colors = [expt.isbuildup[i] ? :orange : :steelblue for i in 1:nslices(expt)]

    ax = Axis(panel[1, 1];
              xlabel="Experiment",
              ylabel="Relative amplitude",
              xticks=(1:nslices(expt), expt.specdata.zlabels))

    hlines!(ax, [0]; linewidth=0)
    hlines!(ax, fit; linewidth=2, color=:red, linestyle=:dash)
    errorbars!(ax, err; whiskerwidth=10)
    return barplot!(ax, x, y; color=colors)
end

function makepeakplot!(gui, state, expt, ::CCRVisualisation)
    @debug "making peak plot for CCR visualisation"

    # Color bars: blue for decay, orange for buildup
    colors = [expt.isbuildup[i] ? :orange : :steelblue for i in 1:nslices(expt)]

    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel="Experiment",
                                 ylabel="Relative amplitude",
                                 xticks=(1:nslices(expt), expt.specdata.zlabels))

    hlines!(ax, [0]; linewidth=0)
    hlines!(ax, state[:peak_plot_fit]; linewidth=2, color=:red, linestyle=:dash)
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    return barplot!(ax, state[:peak_plot_x], state[:peak_plot_y]; color=colors)
end