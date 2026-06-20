"""
    cest2d(inputfilename; B1, Tsat)

Start interactive GUI for analysing 2D CEST (Chemical Exchange Saturation Transfer) data.

# Arguments
- `inputfilename`: NMR data file as a processed data directory containing pseudo-3D data
                   where the first plane is the reference spectrum and subsequent planes
                   are the saturation spectra
- `B1`: Saturation power in Hz
- `Tsat`: Saturation time in seconds

# Example:
```julia
cest2d("path/to/expno/pdata/1"; B1=15, Tsat=0.3)
```
"""
function cest2d(inputfilename; B1, Tsat)
    expt = CESTExperiment(inputfilename, B1, Tsat)
    return gui!(expt)
end

"""
    CESTExperiment <: FixedPeakExperiment

Chemical Exchange Saturation Transfer experiment with reference and saturation spectra.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `frequencies`: Vector of saturation frequencies in ppm
"""
struct CESTExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    frequencies::Vector{Float64}
    B1::Float64
    Tsat::Float64

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function CESTExperiment(specdata, peaks, frequencies, B1, Tsat)
        expt = new(specdata, peaks, frequencies, B1, Tsat,
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

struct CESTVisualisation <: VisualisationStrategy end
visualisationtype(::CESTExperiment) = CESTVisualisation()
primaryparam(::CESTExperiment) = :R1

"""
    CESTExperiment(inputfilename)

Create CEST experiment from a pseudo-3D input file where the first plane is the reference
and the remaining planes are saturation spectra at different frequencies.
"""
function CESTExperiment(inputfilename, B1, Tsat)
    @debug "Creating CEST experiment from $inputfilename with B1=$B1 Hz and Tsat=$Tsat s"
    spec = loadnmr(inputfilename)

    # Extract saturation frequencies from fq3list using proper NMRTools methods
    frequencies = if haskey(acqus(spec), :fq3list)
        fq_list = acqus(spec, :fq3list)
        # Get frequencies in ppm
        ppm(fq_list, dims(spec, F2Dim))
    else
        # Fallback if fq3list is not available
        collect(range(-10.0, 10.0; length=ndims(spec, 3)))
    end

    # Prepare specdata
    specdata = preparespecdata(inputfilename, frequencies, CESTExperiment)
    peaks = Observable(Vector{Peak}())

    return CESTExperiment(specdata, peaks, frequencies, B1, Tsat)
end

# Load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, frequencies, ::Type{CESTExperiment})
    @debug "Preparing spec data for CEST experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)

    # Get 3D data and normalize by scale
    raw_data = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    # Extract slices from the 3D data
    z = eachslice(raw_data; dims=3)

    # Create labels for each saturation frequency
    zlabels = ["Reference"]
    for freq in frequencies[2:end]  # Skip first frequency for reference
        push!(zlabels, "$(round(freq, digits=2)) ppm")
    end

    return SpecData(SingleElementVector(spec),
                    SingleElementVector(x),
                    SingleElementVector(y),
                    z ./ σ,
                    SingleElementVector(1),
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::CESTExperiment, initialposition::Point2f, label="",
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

    # Add post-parameters for CEST analysis
    newpeak.postparameters[:R1] = Parameter("R1", 0.0)
    newpeak.postparameters[:R2] = Parameter("R2", 0.0)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::CESTExperiment, xbounds=nothing, ybounds=nothing)
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

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::CESTExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog = 10

    δsat = expt.frequencies[2:end]
    δ0 = peak.parameters[:y].value[][1]
    R20 = peak.parameters[:R2y].value[][1]
    v0 = 1e-6 * δ0 * expt.specdata.nmrdata[1][2, :bf]
    vsat = 1e-6 * δsat .* expt.specdata.nmrdata[1][2, :bf]
    Tsat = expt.Tsat
    v1 = expt.B1

    zspecobs = peak.parameters[:amp].value[][2:end] ./ peak.parameters[:amp].value[][1]
    p0 = [1.5, R20] # R1, R2 

    # Fit the model
    model(x, p) = Z(x, v0, v1, Tsat, p[1], p[2])
    fit = curve_fit(model, vsat, zspecobs, p0)
    pfit = coef(fit)
    perr = stderror(fit)
    # pfit = p0
    # perr = [0.1, 0.1]

    # Update post-parameters with fitted values
    peak.postparameters[:R1].value[] .= pfit[1]
    peak.postparameters[:R1].uncertainty[] .= perr[1]
    peak.postparameters[:R2].value[] .= pfit[2]
    peak.postparameters[:R2].uncertainty[] .= perr[2]

    @debug "Fitted parameters: $(peak.postparameters)"
    return peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::CESTExperiment, idx)
    if idx == 1
        "Reference"
    else
        "Saturation at $(round(expt.frequencies[idx], digits=2)) ppm ($idx of $(nslices(expt)))"
    end
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::CESTExperiment, idx)
    if idx == 0
        return "No peak selected"
    end

    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "R1: $(peak.postparameters[:R1].value[][1] ± peak.postparameters[:R1].uncertainty[][1]) s⁻¹\n" *
               "R2: $(peak.postparameters[:R2].value[][1] ± peak.postparameters[:R2].uncertainty[][1]) s⁻¹\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::CESTExperiment)
    return "Analysis type: CEST (Chemical Exchange Saturation Transfer)\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Number of planes: $(nslices(expt))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

## Visualisation
function get_cest_data(peak, expt::CESTExperiment)
    @debug "getting CEST data"
    isnothing(peak) && return (Point2f[], [(0.0, 0.0, 0.0)], 0.0, Point2f[])

    # X-axis will be frequency values
    x = expt.frequencies[2:end]

    # Get amplitudes and reference amplitude
    amp = peak.parameters[:amp].value[]
    amp_err = peak.parameters[:amp].uncertainty[]
    ref_amp = 1.0  # Reference amplitude is always 1.0
    ref_err = amp[1] > 0 ? amp_err[1] / amp[1] : 0.1

    # Calculate relative intensities
    y = amp[2:end] ./ amp[1]
    yerr = amp_err[2:end] ./ amp[1]
    # if error > 100% set to 100%
    yerr = min.(yerr, 1.0)

    # Create error tuples for plotting (without propagating reference uncertainty)
    obs_points = Point2f.(x, y)
    obs_err = collect(zip(x, y, yerr))

    # Calculate fit line if peak has been fitted
    if peak.postfitted[]
        δ0 = peak.parameters[:y].value[][1]
        v0 = 1e-6 * δ0 * expt.specdata.nmrdata[1][2, :bf]
        vsat = 1e-6 * x * expt.specdata.nmrdata[1][2, :bf]
        Tsat = expt.Tsat
        v1 = expt.B1
        R1 = peak.postparameters[:R1].value[][1]
        R2 = peak.postparameters[:R2].value[][1]

        ypred = Z.(vsat, v0, v1, Tsat, R1, R2)
        fit_points = Point2f.(x, ypred)
    else
        fit_points = Point2f[]
    end

    return (obs_points, obs_err, ref_err, fit_points)
end

function completestate!(state, expt, ::CESTVisualisation)
    @debug "completing state for CEST visualisation"
    state[:peak_plot_data] = lift(peak -> get_cest_data(peak, expt), state[:current_peak])
    # state[:peak_plot_x] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_obs] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_err] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_ref_err] = lift(d -> d[3], state[:peak_plot_data])
    return state[:peak_plot_fit] = lift(d -> d[4], state[:peak_plot_data])
end

function plot_peak!(panel, peak, expt, ::CESTVisualisation)
    @debug "plotting peak for CEST visualisation"

    obs_points, obs_err, ref_err, fit_points = get_cest_data(peak, expt)

    ax = Axis(panel[1, 1];
              xlabel="Saturation frequency (ppm)",
              ylabel="Relative intensity (I/I₀)")

    hlines!(ax, [0.0]; linewidth=0, color=:black) # sneaky way to ensure axis goes to zero
    hlines!(ax, [1.0]; linewidth=1, color=:gray, linestyle=:dash) # Add reference line at y=1
    hspan!(ax, 1 - ref_err, 1 + ref_err; color=(:gray, 0.2))
    errorbars!(ax, obs_err; whiskerwidth=10)
    scatter!(ax, obs_points)
    return lines!(ax, fit_points; color=:red)
end

function makepeakplot!(gui, state, expt, ::CESTVisualisation)
    @debug "making peak plot for CEST visualisation"
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel="Saturation frequency (ppm)",
                                 ylabel="Relative intensity (I/I₀)",
                                 title="CEST Profile")

    hlines!(ax, [0.0]; linewidth=0, color=:black) # sneaky way to ensure axis goes to zero
    hlines!(ax, [1.0]; linewidth=1, color=:gray, linestyle=:dash)
    hspan!(ax, lift(e -> 1 - e, state[:peak_plot_ref_err]),
           lift(e -> 1 + e, state[:peak_plot_ref_err]);
           color=(:gray, 0.2))
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    scatter!(ax, state[:peak_plot_obs])
    return lines!(ax, state[:peak_plot_fit]; color=:red)
end

## CEST fitting (exchange free)
function Z(vsat, v0, v1, Tsat, R1, R2)
    # https://iopscience.iop.org/article/10.1088/0031-9155/58/22/R221

    ωsat = 2π * vsat # s-1
    ω0 = 2π * v0 # s-1
    ω1 = 2π * v1 # s-1
    Δω = ω0 .- ωsat

    cosθ = @. Δω / sqrt(ω1^2 + Δω^2)
    R1res = cosθ * R1
    Pz = cosθ
    Pzeff = cosθ

    Reff = @. (R2 - R1) * ω1^2 / (ω1^2 + Δω^2)

    Zss = @. Pz * R1res / Reff
    Z = @. (Zss + (Pzeff * Pz - Zss) * exp(-Tsat * Reff)) * exp(-Tsat * R1)

    return Z
end