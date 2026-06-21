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

"""
    trackmaximum(expt, i, xc, yc) -> (x, y)

Position of the most intense point in plane `i` within the peak radius of `(xc, yc)`. Used to
follow a peak from plane to plane during tracking; returns `(xc, yc)` unchanged if the search
window falls off the spectrum.
"""
function trackmaximum(expt::MovingPeakExperiment, i, xc, yc)
    x = expt.specdata.x[i]
    y = expt.specdata.y[i]
    z = expt.specdata.z[i]
    xi = findall(xc - expt.xradius[] .≤ x .≤ xc + expt.xradius[])
    yi = findall(yc - expt.yradius[] .≤ y .≤ yc + expt.yradius[])
    (isempty(xi) || isempty(yi)) && return (xc, yc)
    idx = argmax(@view z[xi, yi])
    return (x[xi[idx[1]]], y[yi[idx[2]]])
end

# Set both the fitted and initial position of one plane (so it displays and seeds the refit).
function setpeakposition!(peak::Peak, i, x, y)
    peak.parameters[:x].initialvalue[][i] = x
    peak.parameters[:x].value[][i] = x
    peak.parameters[:y].initialvalue[][i] = y
    peak.parameters[:y].value[][i] = y
    return
end

"""
    addandtrackpeak!(expt, initialposition, label="")

Add a peak at `initialposition` and track it across every plane by following the local maximum:
the current plane is anchored at the click, then the position is propagated outward in both
directions, each plane seeded from its neighbour's tracked position. The subsequent fit then
refines each plane within its radius. Good for non-crowded series (titrations); for crowded
regions, add with `A` and adjust planes by hand instead.
"""
function addandtrackpeak!(expt::MovingPeakExperiment, initialposition, label="")
    addpeak!(expt, initialposition, label)
    peak = expt.peaks[][end]
    s = expt.state[][:current_slice][]

    # Anchor at the current plane (refine the click to the local maximum), then propagate.
    xc, yc = trackmaximum(expt, s, initialposition[1], initialposition[2])
    setpeakposition!(peak, s, xc, yc)
    px, py = xc, yc
    for i in (s + 1):nslices(expt)
        px, py = trackmaximum(expt, i, px, py)
        setpeakposition!(peak, i, px, py)
    end
    px, py = xc, yc
    for i in (s - 1):-1:1
        px, py = trackmaximum(expt, i, px, py)
        setpeakposition!(peak, i, px, py)
    end

    peak.touched[] = true
    notify(expt.peaks)
    return length(expt.peaks[])
end

"""Simulate single peak according to its per-plane position, linewidth and amplitude."""
function simulate!(z, peak::Peak, expt::MovingExperiment, xbounds=nothing, ybounds=nothing)
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
        # Scale by this plane's OWN R2x/R2y so the fitted amplitude is the peak height,
        # decoupled from linewidth. Unlike fixed-peak experiments (where linewidths are shared),
        # each plane must use its own linewidths - otherwise plane 1's linewidth scales every
        # plane's amplitude and couples the whole peak's fit together.
        zx = NMRTools.NMRBase._lineshape(2π * hz(x0, xaxis), R2x, 2π * hz(xs, xaxis),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x * R2y) *
             NMRTools.NMRBase._lineshape(2π * hz(y0, yaxis), R2y, 2π * hz(ys, yaxis),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

# Postfit: NoFitting just marks the peak fitted (the generic Experiment fallback in
# experiments.jl already does this); position-based physical models are added separately.

function addpeakhint(::MovingPeakExperiment)
    return "Press (A) to add a peak, or (T) to add and track it across planes, under the cursor"
end

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
    n = nslices(expt)
    length(expt.specdata.zlabels) == 1 && return "Slice $idx of $n"
    lbl = string(expt.specdata.zlabels[idx])
    length(lbl) > 20 && (lbl = first(lbl, 20) * "…")
    return "$lbl ($idx of $n)"
end

"""
    add_moving_overlays!(g, state, expt::MovingPeakExperiment)

Add the moving-peak overlays to the contour panel:
- a faint polyline tracing each peak's fitted position across all planes (the trajectory),
  so the walk is visible at a glance;
- a toggleable "Context" overlay drawing every plane's contours faintly behind the current
  one, for orienting peaks that move a long way.
"""
function add_moving_overlays!(g, state, expt::MovingPeakExperiment)
    ax = g[:axcontour]

    # --- per-peak position trajectories across all planes ---
    # One flat point list with NaN separators between peaks, so a single lines! call draws
    # every trajectory as disconnected segments, plus matching per-vertex dots. Colours mirror
    # the peak markers: touched/unfitted → red, fitted → blue, selected → lime. Concrete RGBAf
    # values are used (a Vector{Symbol} colour is not honoured by lines!).
    statuscolour(j, sel, touched) = j == sel ? RGBAf(0, 1, 0, 1) :
                                    touched ? RGBAf(1, 0, 0, 1) : RGBAf(0, 0, 1, 1)
    # Single source of truth, recomputed whenever the peaks (positions/fit status) or the
    # selection change; the points and colours derive from it so they stay length-consistent.
    state[:trajectorydata] = lift(expt.peaks, state[:current_peak_idx]) do peaks, sel
        pts = Point2f[]
        cols = RGBAf[]
        for (j, peak) in enumerate(peaks)
            c = statuscolour(j, sel, peak.touched[])
            xs = peak.parameters[:x].value[]
            ys = peak.parameters[:y].value[]
            for i in 1:length(xs)
                push!(pts, Point2f(xs[i], ys[i]))
                push!(cols, c)
            end
            push!(pts, Point2f(NaN32, NaN32))  # break between peaks
            push!(cols, c)
        end
        return (pts, cols)
    end
    # Register the colour derivation before the position one so it is up to date by the time a
    # position change triggers a redraw.
    state[:trajectorycolours] = lift(d -> d[2], state[:trajectorydata])
    state[:trajectories] = lift(d -> d[1], state[:trajectorydata])

    g[:plttrajectories] = lines!(ax, state[:trajectories];
                                 color=state[:trajectorycolours], linewidth=2.5)
    g[:plttrajectorypts] = scatter!(ax, state[:trajectories];
                                    color=state[:trajectorycolours], markersize=7)

    # --- faint context: every plane's contours, off by default ---
    # Drawn above the (opaque white) mask heatmap so they remain visible; the small positive z
    # keeps them under the peak markers and trajectory. The "Show all" toggle widget is created
    # in gui! (just left of the Fitting toggle); here we attach its plots and handler.
    g[:pltotherplanes] = map(1:nslices(expt)) do i
        p = contour!(ax, expt.specdata.x[i], expt.specdata.y[i], expt.specdata.z[i];
                     levels=g[:contourlevels], color=(:grey60, 0.3), visible=false)
        translate!(p, 0, 0, 1)
        return p
    end
    on(g[:toggleother].active) do active
        for p in g[:pltotherplanes]
            p.visible[] = active
        end
    end

    return nothing
end

# Legend label for a plane: its independent-variable value when one was supplied,
# otherwise the plane index.
function _planelabel(expt::MovingExperiment, i)
    expt.x == collect(1.0:nslices(expt)) && return "Plane $i"
    return string(round(expt.x[i]; digits=3))
end

"""
    summaryplot(expt::MovingExperiment; weights=(1.0, 0.14), title, size, include_unassigned)

Default summary for a moving-peak experiment: combined chemical-shift perturbation Δδ
against residue number, with one series per plane beyond the first (each plane's shift is
measured relative to plane 1).

Δδ = √((w₁·Δδ_x)² + (w₂·Δδ_y)²), with `weights` weighting the two dimensions. The default
`(1.0, 0.14)` assumes the F1 (x) axis is ¹H and the F2 (y) axis is ¹⁵N (Williamson 2013); pass
`weights` to suit other nuclei (e.g. `(1.0, 0.25)` for ¹³C).
"""
function summaryplot(expt::MovingExperiment; weights=(1.0, 0.14), title="", size=nothing,
                     include_unassigned=false)
    n = nslices(expt)
    peaks = sortedpeaks(expt)
    if !include_unassigned
        assigned = filter(p -> extract_residue_number(p.label[]) > 0, peaks)
        isempty(assigned) || (peaks = assigned)
    end

    figkw = isnothing(size) ? NamedTuple() : (; size=size)
    fig = Figure(; figkw...)
    ax = Axis(fig[1, 1]; xlabel="Residue number", ylabel="Δδ / ppm", title=title,
              xgridvisible=false, ygridvisible=false)

    resnums = Float64[extract_residue_number(p.label[]) for p in peaks]
    wx, wy = weights
    for i in 2:n
        Δδ = map(peaks) do peak
            dx = peak.parameters[:x].value[][i] - peak.parameters[:x].value[][1]
            dy = peak.parameters[:y].value[][i] - peak.parameters[:y].value[][1]
            return sqrt((wx * dx)^2 + (wy * dy)^2)
        end
        scatterlines!(ax, resnums, Float64.(Δδ); label=_planelabel(expt, i))
    end
    hlines!(ax, [0]; linewidth=0)  # invisible: forces zero into the y-range
    n > 2 && axislegend(ax; position=:lt, framevisible=false)

    return fig
end
