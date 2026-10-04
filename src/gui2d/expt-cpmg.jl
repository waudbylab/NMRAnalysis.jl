"""
    cpmg2d(inputfilename; Trelax, vCPMG, skipplanes=nothing, prompt=isinteractive())
    cpmg2d(inputfilename; Trelax, ncyc, skipplanes=nothing, prompt=isinteractive())
        -> CPMGExperiment

Start interactive GUI for analysing 2D CPMG relaxation dispersion data. The window blocks
until it is closed, and the analysis is returned (see [`results`](@ref)).

The CPMG frequencies are taken from `vCPMG`, else from `ncyc`, else from the cycle numbers
in the `vclist`; `Trelax` is taken from its keyword. Anything missing is asked for.

# Arguments
- `inputfilename`: NMR data file as a processed data directory containing pseudo-3D data
                   where the first plane is the reference spectrum and subsequent planes
                   are the saturation spectra
- `Trelax`: Relaxation time in seconds
- `vCPMG`: list of CPMG frequencies in Hz, use zero for reference spectrum
- `ncyc`: list of CPMG cycle numbers, use zero for reference spectrum.
          When provided, vCPMG is calculated as ncyc/Trelax
- `skipplanes`: Planes (1-based) left out of the fitting. At least one reference plane must
  remain.

# Examples:
```julia
# Direct specification of CPMG frequencies
cpmg2d("path/to/expno"; Trelax=0.04, vCPMG=[0, 25, 50, 75, 100])

# Using cycle numbers (vCPMG calculated automatically)
ncyc = [0, 1, 2, 3, 4]
cpmg2d("path/to/expno"; Trelax=0.04, ncyc=ncyc)

# Cycle numbers from the vclist
cpmg2d("path/to/expno"; Trelax=0.04)
```
"""
function cpmg2d(inputfilename; Trelax=nothing, vCPMG=nothing, ncyc=nothing,
                skipplanes=nothing, peaklist=nothing, prompt::Bool=isinteractive())
    isnothing(vCPMG) || isnothing(ncyc) ||
        throw(ArgumentError("Cannot specify both vCPMG and ncyc"))
    given = inputfilename
    inputfilename = asexptpath(inputfilename)
    spec = loadnmr(inputfilename)
    n = size(spec, 3)
    Trelax = @something(Trelax, ask("CPMG relaxation time Trelax"; unit="s", prompt))
    vCPMG = @something(isnothing(vCPMG) ? nothing : readvalues(vCPMG),
                       isnothing(ncyc) ? nothing : readvalues(ncyc) ./ Trelax,
                       let c = acqusvalue(spec, :vclist)
                           isnothing(c) ? nothing : collect(Float64, c) ./ Trelax
                       end,
                       askvector("CPMG cycle numbers", n; prompt) ./ Trelax)
    expt = CPMGExperiment(inputfilename, Trelax, vCPMG;
                          skipplanes=checkskipplanes(skipplanes, n))
    call = analysiscall("cpmg2d", given; Trelax, vCPMG, skipplanes=nonempty(expt.skipplanes))
    return gui!(expt; peaklist, call)
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
    skipplanes::Vector{Int}

    function CPMGExperiment(specdata, peaks, Trelax, vCPMG; skipplanes=Int[])
        any(i -> vCPMG[i] ≈ 0 && i ∉ skipplanes, eachindex(vCPMG)) ||
            throw(ArgumentError("no reference plane (vCPMG = 0) is left unskipped"))
        expt = new(specdata, peaks, Trelax, vCPMG,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}(),
                   collect(Int, skipplanes))
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

struct CPMGVisualisation <: VisualisationStrategy end
visualisationtype(::CPMGExperiment) = CPMGVisualisation()
primaryparam(::CPMGExperiment) = :R20

"""
    CPMGExperiment(inputfilename, Trelax, vCPMG; skipplanes=Int[])

Create CPMG experiment from a pseudo-3D input file. Zero frequency in `vCPMG` indicates the reference plane.
"""
function CPMGExperiment(inputfilename, Trelax, vCPMG; skipplanes=Int[])
    @debug "Creating CPMG experiment from $inputfilename with Trelax=$Trelax s and vCPMG=$vCPMG Hz"
    spec = loadnmr(inputfilename)

    # check that vCPMG is a vector of frequencies with length matching spectrum size
    if length(vCPMG) != size(spec, X3Dim)
        error("vCPMG must be a vector of frequencies with length matching the number of slices in the spectrum.")
    end

    # Prepare specdata
    specdata = preparespecdata(inputfilename, vCPMG, CPMGExperiment)
    peaks = Observable(Vector{Peak}())

    return CPMGExperiment(specdata, peaks, Trelax, vCPMG; skipplanes)
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

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::CPMGExperiment)
    @debug "Post-fitting peak $(peak.label)" #maxlog = 10

    vCPMG = expt.vCPMG
    Trelax = expt.Trelax
    used = [i ∉ skipset(expt) for i in eachindex(vCPMG)]
    refidx = (vCPMG .≈ 0) .& used
    cpmgidx = .!(vCPMG .≈ 0) .& used
    cpmglist = vCPMG[cpmgidx]
    R20 = peak.parameters[:R2y].value[][1]

    Aref = peak.parameters[:amp].value[][refidx] .±
           peak.parameters[:amp].uncertainty[][refidx]
    Aref = sum(Aref) / length(Aref) # average reference amplitude

    Acpmg = peak.parameters[:amp].value[][cpmgidx] .±
            peak.parameters[:amp].uncertainty[][cpmgidx]
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
    perr = stderrors(fit)
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
    used = [i ∉ skipset(expt) for i in eachindex(vCPMG)]
    refidx = (vCPMG .≈ 0) .& used
    cpmgidx = .!(vCPMG .≈ 0) .& used
    x = vCPMG[cpmgidx]

    Aref = peak.parameters[:amp].value[][refidx] .±
           peak.parameters[:amp].uncertainty[][refidx]
    Aref = sum(Aref) / length(Aref) # average reference amplitude

    Acpmg = peak.parameters[:amp].value[][cpmgidx] .±
            peak.parameters[:amp].uncertainty[][cpmgidx]
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
              ylabel="R₂,eff (s⁻¹)")

    hlines!(ax, [0.0]; color=:grey50, linewidth=1)
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

    hlines!(ax, [0.0]; color=:grey50, linewidth=1)
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    scatter!(ax, state[:peak_plot_obs])
    return lines!(ax, state[:peak_plot_fit]; color=:red)
end
