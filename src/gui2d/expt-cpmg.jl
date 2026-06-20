"""
    cpmg2d(inputfilename; Trelax, vCPMG)
    cpmg2d(inputfilename; Trelax, ncyc)

Start interactive GUI for analysing 2D CPMG relaxation dispersion data.

# Arguments
- `inputfilename`: NMR data file as a processed data directory containing pseudo-3D data
                   where the first plane is the reference spectrum and subsequent planes
                   are the saturation spectra
- `Trelax`: Relaxation time in seconds
- `vCPMG`: list of CPMG frequencies in Hz, use zero for reference spectrum
- `ncyc`: list of CPMG cycle numbers, use zero for reference spectrum. 
          When provided, vCPMG is calculated as ncyc/Trelax

# Examples:
```julia
# Direct specification of CPMG frequencies
cpmg2d("path/to/expno"; Trelax=0.04, vCPMG=[0, 25, 50, 75, 100])

# Using cycle numbers (vCPMG calculated automatically)
ncyc = [0, 1, 2, 3, 4]
cpmg2d("path/to/expno"; Trelax=0.04, ncyc=ncyc)
```
"""
function cpmg2d(inputfilename; Trelax, vCPMG=nothing, ncyc=nothing)
    if !isnothing(vCPMG) && !isnothing(ncyc)
        throw(ArgumentError("Cannot specify both vCPMG and ncyc"))
    elseif !isnothing(ncyc)
        vCPMG = ncyc ./ Trelax
    elseif isnothing(vCPMG)
        throw(ArgumentError("Must specify either vCPMG or ncyc"))
    end

    expt = CPMGExperiment(inputfilename, Trelax, vCPMG)
    return gui!(expt)
end

"""
    CPMGExperiment <: FixedPeakExperiment

CPMG experiment with reference plane and relaxation planes.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `vCPMG`: Vector of CPMG frequencies in Hz (zero for reference)
- `Trelax`: Relaxation time in seconds
"""
struct CPMGExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    Trelax::Float64
    vCPMG::Vector{Float64}

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function CPMGExperiment(specdata, peaks, Trelax, vCPMG)
        expt = new(specdata, peaks, Trelax, vCPMG,
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

struct CPMGVisualisation <: VisualisationStrategy end
visualisationtype(::CPMGExperiment) = CPMGVisualisation()
primaryparam(::CPMGExperiment) = :R20

"""
    CPMGExperiment(inputfilename, Trelax, vCPMG)

Create CPMG experiment from a pseudo-3D input file. Zero frequency in `vCPMG` indicates the reference plane.
"""
function CPMGExperiment(inputfilename, Trelax, vCPMG)
    @debug "Creating CPMG experiment from $inputfilename with Trelax=$Trelax s and vCPMG=$vCPMG Hz"
    spec = loadnmr(inputfilename)

    # check that vCPMG is a vector of frequencies with length matching spectrum size
    if length(vCPMG) != size(spec, X3Dim)
        error("vCPMG must be a vector of frequencies with length matching the number of slices in the spectrum.")
    end

    # Prepare specdata
    specdata = preparespecdata(inputfilename, vCPMG, CPMGExperiment)
    peaks = Observable(Vector{Peak}())

    return CPMGExperiment(specdata, peaks, Trelax, vCPMG)
end

# Load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, vCPMG, ::Type{CPMGExperiment})
    @debug "Preparing spec data for CPMG experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)

    # Get 3D data and normalize by scale
    raw_data = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    # Extract slices from the 3D data
    z = eachslice(raw_data; dims=3)

    # Create labels for each saturation frequency
    zlabels = ["$(round(v, digits=0)) Hz" for v in vCPMG]

    return SpecData(SingleElementVector(spec),
                    SingleElementVector(x),
                    SingleElementVector(y),
                    z ./ σ,
                    SingleElementVector(1),
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::CPMGExperiment, initialposition::Point2f, label="",
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

    # Add post-parameters for CPMG analysis
    newpeak.postparameters[:R20] = Parameter("R2,0", 10.0)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::CPMGExperiment, xbounds=nothing, ybounds=nothing)
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
function postfit!(peak::Peak, expt::CPMGExperiment)
    @debug "Post-fitting peak $(peak.label)" #maxlog = 10

    vCPMG = expt.vCPMG
    Trelax = expt.Trelax
    refidx = vCPMG .≈ 0
    cpmglist = vCPMG[.!refidx]
    R20 = peak.parameters[:R2y].value[][1]

    Aref = peak.parameters[:amp].value[][refidx] .±
           peak.parameters[:amp].uncertainty[][refidx]
    Aref = sum(Aref) / length(Aref) # average reference amplitude

    Acpmg = peak.parameters[:amp].value[][.!refidx] .±
            peak.parameters[:amp].uncertainty[][.!refidx]
    R2effobs = @. log(Acpmg / Aref) / (-Trelax) # effective R2 from CPMG amplitudes
    y = Measurements.value.(R2effobs) # convert to Float64 for fitting
    yerr = Measurements.uncertainty.(R2effobs) # uncertainties in R2eff
    w = 1 ./ yerr .^ 2 # weights for fitting
    p0 = [R20]

    # Fit the model
    @debug "Fitting model to CPMG data" cpmglist R2effobs p0
    # model(x, p) = Z(x, v0, v1, Tsat, p[1], p[2])
    model(x, p) = ones(length(x)) * p[1] # no exchange, R2eff = R20
    fit = curve_fit(model, cpmglist, y, w, p0)
    pfit = coef(fit)
    perr = stderror(fit)
    # pfit = p0
    # perr = [0.1]
    @debug "Fitted parameters: $(pfit), uncertainties: $(perr)"

    # Update post-parameters with fitted values
    peak.postparameters[:R20].value[] .= pfit[1]
    peak.postparameters[:R20].uncertainty[] .= perr[1]

    @debug "Fitted parameters: $(peak.postparameters)"
    return peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::CPMGExperiment, idx)
    if expt.vCPMG[idx] ≈ 0
        "Reference"
    else
        "$(expt.vCPMG[idx]) Hz ($idx of $(nslices(expt)))"
    end
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::CPMGExperiment, idx)
    if idx == 0
        return "No peak selected"
    end

    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "R2,0: $(peak.postparameters[:R20].value[][1] ± peak.postparameters[:R20].uncertainty[][1]) s⁻¹\n" *
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
function experimentinfo(expt::CPMGExperiment)
    return "Analysis type: CPMG (relaxation dispersion)\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Number of planes: $(nslices(expt))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

## Visualisation
function get_cpmg_data(peak, expt::CPMGExperiment)
    @debug "getting CPMG data"
    isnothing(peak) && return (Point2f[], [(0.0, 0.0, 0.0)], Point2f[])

    vCPMG = expt.vCPMG
    Trelax = expt.Trelax
    refidx = vCPMG .≈ 0
    x = vCPMG[.!refidx]
    # @debug "CPMG frequencies (Hz)" x

    Aref = peak.parameters[:amp].value[][refidx] .±
           peak.parameters[:amp].uncertainty[][refidx]
    Aref = sum(Aref) / length(Aref) # average reference amplitude

    Acpmg = peak.parameters[:amp].value[][.!refidx] .±
            peak.parameters[:amp].uncertainty[][.!refidx]
    R2effobs = @. log(Acpmg / Aref) / (-Trelax) # effective R2 from CPMG amplitudes
    y = Measurements.value.(R2effobs) # convert to Float64 for plotting
    yerr = Measurements.uncertainty.(R2effobs) # uncertainties in R2eff

    # Create error tuples for plotting (without propagating reference uncertainty)
    obs_points = Point2f.(x, y)
    obs_err = collect(zip(x, y, yerr))

    # Calculate fit line if peak has been fitted
    if peak.postfitted[]
        ypred = peak.postparameters[:R20].value[][1]
        fit_points = [Point2f(xval, ypred) for xval in sort(x)]
    else
        fit_points = Point2f[]
    end
    @debug "get_cpmg_data" obs_points obs_err fit_points

    return (obs_points, obs_err, fit_points)
end

function completestate!(state, expt, ::CPMGVisualisation)
    @debug "completing state for CPMG visualisation"
    state[:peak_plot_data] = lift(peak -> get_cpmg_data(peak, expt),
                                  state[:current_peak])
    # # state[:peak_plot_x] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_obs] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_err] = lift(d -> d[2], state[:peak_plot_data])
    return state[:peak_plot_fit] = lift(d -> d[3], state[:peak_plot_data])
end

function plot_peak!(panel, peak, expt, ::CPMGVisualisation)
    @debug "plotting peak for CPMG visualisation"

    obs_points, obs_err, fit_points = get_cpmg_data(peak, expt)

    ax = Axis(panel[1, 1];
              xlabel="νCPMG (Hz)",
              ylabel="R₂,eff (s⁻¹)",)

    hlines!(ax, [0.0]; linewidth=0, color=:black) # sneaky way to ensure axis goes to zero
    errorbars!(ax, obs_err; whiskerwidth=10)
    scatter!(ax, obs_points)
    return lines!(ax, fit_points; color=:red)
end

function makepeakplot!(gui, state, expt, ::CPMGVisualisation)
    @debug "making peak plot for CPMG visualisation"
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel="νCPMG (Hz)",
                                 ylabel="R₂,eff (s⁻¹)",
                                 title="CPMG Profile")

    hlines!(ax, [0.0]; linewidth=0, color=:black) # sneaky way to ensure axis goes to zero
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    scatter!(ax, state[:peak_plot_obs])
    return lines!(ax, state[:peak_plot_fit]; color=:red)
end
