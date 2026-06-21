"""
    MovingExperiment <: MovingPeakExperiment

Generic experiment for measurement of peak *positions* that change across a series of 2D
spectra. Unlike [`IntensityExperiment`](@ref), where each peak sits at one chemical shift and
only its amplitude is read out, here each peak's position (and linewidth) is an independent
fit parameter in every plane. This is the common core for titrations, coupling-constant and
RDC measurements: the spectrum fit recovers a per-plane position trajectory, and a postfit
`model` turns that trajectory into a physical quantity.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks (each holds per-plane `:x`, `:y`, `:R2x`, `:R2y`, `:amp`)
- `x`: Per-plane independent variable (e.g. ligand concentration); defaults to plane index
- `model`: Postfit model relating the position trajectory to the independent variable
- `visualisation`: Strategy for the per-peak plot panel
"""
struct MovingExperiment <: MovingPeakExperiment
    specdata::Any
    peaks::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    x::Vector{Float64}
    model::FittingModel
    visualisation::VisualisationStrategy

    function MovingExperiment(specdata, peaks, model=NoFitting(), xvalues=nothing,
                              visualisation=CrossSectionVisualisation())
        if isnothing(xvalues)
            xvalues = 1.0 * collect(1:length(specdata.z))
        end
        expt = new(specdata, peaks,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}(),
                   xvalues, model, visualisation)
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

visualisationtype(expt::MovingExperiment) = expt.visualisation

"""
    movingfit2d(inputfilenames, xvalues=nothing)

Start an interactive GUI for analysing a series of 2D spectra in which **peak positions
change** from plane to plane (e.g. a titration, or a coupling-constant / RDC measurement).
Each peak is fitted to a 2D Lorentzian lineshape independently in every plane, so its
position and linewidth are free to move; no physical model is applied to the trajectory.

Use this to measure how peak positions and linewidths evolve across a series for downstream
analysis, or as the basis for the position-based physical models (titration, coupling).

# Arguments
- `inputfilenames`: A single path string (pseudo-3D dataset) or a vector of path strings
  (one file per plane) pointing to processed Bruker data directories.
- `xvalues`: Optional vector of independent-variable values (one per plane), or a string
  giving a path to a text file containing them (one per line; `#` comments ignored). When
  omitted, the plane index is used.

# Example
```julia
# A titration series, with ligand concentrations
movingfit2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"], [0.0, 0.5, 1.0])
```
"""
function movingfit2d(inputfilenames, xvalues=nothing)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    xval = if xvalues isa String
        Float64.(vec(readdlm(xvalues; comments=true)))
    elseif xvalues isa Vector
        Float64.(xvalues)
    else
        nothing
    end

    expt = MovingExperiment(specdata, peaks, NoFitting(), xval)

    return gui!(expt)
end

"""Add peak to experiment, setting up per-plane position, linewidth and amplitude parameters."""
function addpeak!(expt::MovingExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add moving peak $label at $initialposition"

    n = nslices(expt)
    # Seed every plane at the clicked position - a StandardVector position makes :x/:y
    # independent per-plane fit parameters (each bounded within ±radius of its own value).
    positions = fill(initialposition, n)
    newpeak = Peak(positions, label, xradius, yradius)

    # per-plane linewidths and amplitudes
    R2x = Parameter("R2x", MaybeVector(fill(30.0, n)); minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", MaybeVector(fill(15.0, n)); minvalue=1.0, maxvalue=100.0)

    x0, y0 = initialposition
    amp0 = map(1:n) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)

    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    # Post-parameters based on model type (none for NoFitting)
    setup_post_parameters!(newpeak, expt.model)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to its per-plane position, linewidth and amplitude."""
function simulate!(z, peak::Peak, expt::MovingExperiment, xbounds=nothing, ybounds=nothing)
    # Reference linewidths (plane 1) used only to decouple amplitude from linewidth scaling.
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]

    for i in 1:nslices(expt)
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
        zy = (π^2 * amp * R2x0 * R2y0) *
             NMRTools.NMRBase._lineshape(2π * hz(y0, yaxis), R2y, 2π * hz(ys, yaxis),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

# Postfit: NoFitting just marks the peak fitted (the generic Experiment fallback in
# experiments.jl already does this); position-based physical models are added separately.

function peakinfotext(expt::MovingExperiment, idx)
    idx == 0 && return "No peak selected"
    peak = expt.peaks[][idx]
    if !peak.postfitted[]
        return "Peak: $(peak.label[])\nNot fitted"
    end

    n = nslices(expt)
    x = peak.parameters[:x].value[]
    y = peak.parameters[:y].value[]
    info = ["Peak: $(peak.label[])",
            "",
            "Plane 1: δX $(round(x[1]; digits=3)), δY $(round(y[1]; digits=3)) ppm",
            "Plane $n: δX $(round(x[n]; digits=3)), δY $(round(y[n]; digits=3)) ppm",
            "",
            "ΔδX (plane $n − 1): $(round(x[n] - x[1]; digits=4)) ppm",
            "ΔδY (plane $n − 1): $(round(y[n] - y[1]; digits=4)) ppm"]
    return join(info, "\n")
end

function experimentinfo(expt::MovingExperiment)
    info = ["Analysis type: Moving peak",
            "Model: $(typeof(expt.model))",
            "Filename: $(expt.specdata.nmrdata[1][:filename])",
            "Number of peaks: $(length(expt.peaks[]))",
            "Number of planes: $(nslices(expt))"]
    return join(info, "\n")
end

function slicelabel(expt::MovingExperiment, idx)
    if length(expt.specdata.zlabels) == 1
        "Slice $idx of $(nslices(expt))"
    else
        "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
    end
end
