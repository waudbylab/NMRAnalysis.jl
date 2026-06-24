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
    x0, y0 = initialposition
    # Initial linewidths from the radii: take half the radius as a FWHM estimate and convert to
    # an R2 rate (R2 = π·FWHM in Hz), clamped to the fit bounds. Adapts to each dimension's
    # frequency scale rather than using a fixed guess.
    xaxis = dims(expt.specdata.nmrdata[1], F1Dim)
    yaxis = dims(expt.specdata.nmrdata[1], F2Dim)
    R2x = Parameter("R2x", MaybeVector(fill(r2_from_radius(xaxis, xradius, x0), n));
                    minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", MaybeVector(fill(r2_from_radius(yaxis, yradius, y0), n));
                    minvalue=1.0, maxvalue=100.0)

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

# Set both the fitted and initial position of one plane (so it displays and seeds the refit),
# and resample the amplitude at the new position.
function setpeakposition!(expt::MovingPeakExperiment, peak::Peak, i, x, y)
    peak.parameters[:x].initialvalue[][i] = x
    peak.parameters[:x].value[][i] = x
    peak.parameters[:y].initialvalue[][i] = y
    peak.parameters[:y].value[][i] = y
    reinitialise_amplitude!(expt, peak, i)
    return
end

# R2 rate (s⁻¹) for an initial linewidth guess: half the search radius taken as a FWHM, converted
# from Hz (R2 = π·FWHM), clamped to the fit bounds.
function r2_from_radius(axis, radius, pos)
    Δhz = abs(hz(pos + radius, axis) - hz(pos, axis))
    return clamp(π * Δhz / 2, 1.0, 100.0)
end

# Resample plane `i`'s amplitude from the spectrum at the peak's current initial position. After
# the per-plane R2 scaling in simulatepeakplane!, the fitted amplitude is the peak height, so the
# nearest grid intensity is a good seed - kept current whenever the position moves.
function reinitialise_amplitude!(expt::MovingPeakExperiment, peak::Peak, i)
    x0 = peak.parameters[:x].initialvalue[][i]
    y0 = peak.parameters[:y].initialvalue[][i]
    ix = findnearest(expt.specdata.x[i], x0)
    iy = findnearest(expt.specdata.y[i], y0)
    a = expt.specdata.z[i][ix, iy]
    peak.parameters[:amp].initialvalue[][i] = a
    peak.parameters[:amp].value[][i] = a
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
    setpeakposition!(expt, peak, s, xc, yc)
    px, py = xc, yc
    for i in (s + 1):nslices(expt)
        px, py = trackmaximum(expt, i, px, py)
        setpeakposition!(expt, peak, i, px, py)
    end
    px, py = xc, yc
    for i in (s - 1):-1:1
        px, py = trackmaximum(expt, i, px, py)
        setpeakposition!(expt, peak, i, px, py)
    end

    peak.touched[] = true
    notify(expt.peaks)
    return length(expt.peaks[])
end

"""Simulate single peak across all planes (used for the displayed fit spectrum)."""
function simulate!(z, peak::Peak, expt::MovingExperiment, xbounds=nothing, ybounds=nothing)
    for i in 1:nslices(expt)
        xb = isnothing(xbounds) ? nothing : xbounds[i]
        yb = isnothing(ybounds) ? nothing : ybounds[i]
        simulatepeakplane!(z[i], peak, expt, i, xb, yb)
    end
end

"""Simulate a single peak into one plane `i`. `xbounds`/`ybounds`, if given, are the boolean
masks restricting `z` to that plane's fit window."""
function simulatepeakplane!(z, peak::Peak, expt::MovingPeakExperiment, i,
                            xbounds=nothing, ybounds=nothing)
    xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
    yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
    x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds]
    y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds]

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
    # Scale by this plane's OWN R2x/R2y so the fitted amplitude is the peak height, decoupled
    # from linewidth. Unlike fixed-peak experiments (where linewidths are shared), each plane
    # must use its own linewidths - otherwise plane 1's linewidth scales every plane's amplitude.
    zx = NMRTools.NMRBase._lineshape(2π * hz(x0, xaxis), R2x, 2π * hz(xs, xaxis),
                                     xaxis[:window], RealLineshape())
    zy = (π^2 * amp * R2x * R2y) *
         NMRTools.NMRBase._lineshape(2π * hz(y0, yaxis), R2y, 2π * hz(ys, yaxis),
                                     yaxis[:window], RealLineshape())
    z[xi, yi] .+= zx .* zy'
    return z
end

# --- per-plane fitting -------------------------------------------------------
# The per-plane lineshape model is separable across planes, so a moving-peak cluster is fitted
# one plane at a time rather than as one joint optimisation over all 5×nplanes parameters. This
# is both better-conditioned (≈5 params/peak per fit) and fully decoupled: adjusting one plane's
# position cannot disturb another plane's fit. A single joint Levenberg–Marquardt fit couples
# the planes through its shared damping parameter and finite iteration/time budget, so an
# under-converged plane drags the others.

function fit!(cluster::Vector{Int}, expt::MovingPeakExperiment,
              mygen=expt.state[][:fit_generation][], t0=time())
    peaks = [expt.peaks[][i] for i in cluster]
    m = mask(cluster, expt)
    xbounds, ybounds = bounds(m)
    for i in 1:nslices(expt)
        fitplane!(peaks, expt, i, m[i], xbounds[i], ybounds[i], mygen, t0)
    end
    for peak in peaks
        peak.touched.val = false
    end
end

function fitplane!(peaks, expt::MovingPeakExperiment, i, mi, xbi, ybi, mygen, t0)
    p0 = packplane(peaks, i, :initial)
    pmin = packplane(peaks, i, :min)
    pmax = packplane(peaks, i, :max)
    p0 = clamp.(p0, pmin, pmax)

    bigmask = vec(mi)
    smallmask = vec(mi[xbi, ybi])
    zobs = vec(expt.specdata.z[i])[bigmask]
    zbuf = similar(expt.specdata.z[i][xbi, ybi])
    zsimm = similar(zobs)

    function resid(p)
        (expt.state[][:fit_generation][] != mygen) && throw(FitCancelled())
        (time() - t0 > FIT_TIME_BUDGET) && throw(FitCancelled())
        Threads.nthreads() == 1 && yield()
        unpackplane!(copy(p), peaks, i, :value)
        fill!(zbuf, 0.0)
        for peak in peaks
            simulatepeakplane!(zbuf, peak, expt, i, xbi, ybi)
        end
        zsimm .= vec(zbuf)[smallmask]
        return zobs - zsimm
    end

    sol = LsqFit.lmfit(resid, p0, Float64[]; lower=pmin, upper=pmax, autodiff=:finite,
                       maxIter=50, x_tol=1e-4, g_tol=1e-6)
    unpackplane!(coef(sol), peaks, i, :value)
    return unpackplane!(stderror(sol), peaks, i, :uncertainty)
end

# Pack/unpack just plane `i`'s parameters (in OrderedDict order: x, y, R2x, R2y, amp).
function packplane(peaks, i, quantity=:value)
    p = Float64[]
    for peak in peaks
        packplane!(p, peak, i, quantity)
    end
    return p
end

function packplane!(p, peak::Peak, i, quantity)
    for (sym, par) in peak.parameters
        v = if quantity == :initial
            par.initialvalue[][i]
        elseif quantity == :min
            planebound(peak, sym, par, i, :min)
        elseif quantity == :max
            planebound(peak, sym, par, i, :max)
        else
            par.value[][i]
        end
        push!(p, v)
    end
    return p
end

# Position (:x/:y) is bounded within ±radius of this plane's own initial value; other
# parameters keep their fixed scalar bounds (R2: 1–100 s⁻¹; amplitude: unbounded).
function planebound(peak, sym, par, i, which)
    if sym === :x || sym === :y
        r = sym === :x ? peak.xradius[] : peak.yradius[]
        return which === :min ? par.initialvalue[][i] - r : par.initialvalue[][i] + r
    end
    b = which === :min ? par.minvalue[] : par.maxvalue[]
    return b isa AbstractVector ? b[i] : b
end

function unpackplane!(v, peaks, i, quantity=:value)
    for peak in peaks
        unpackplane!(v, peak, i, quantity)
    end
end

function unpackplane!(v, peak::Peak, i, quantity)
    for (_, par) in peak.parameters
        val = popfirst!(v)
        if quantity == :value
            par.value[][i] = val
        elseif quantity == :uncertainty
            par.uncertainty[][i] = val
        end
    end
end

# --- postfit ----------------------------------------------------------------
# Moving-peak postfit dispatches on the model. The default (NoFitting and any model without a
# specific method) just marks the peak fitted; position-based physical models override.
postfit!(peak::Peak, expt::MovingExperiment) = postfitmoving!(peak, expt, expt.model)
function postfitmoving!(peak::Peak, ::MovingExperiment, ::FittingModel)
    peak.postfitted[] = true
    return
end

# --- RDC / coupling measurement ---------------------------------------------

"""
    RDCModel

Derive a scalar coupling `J` and residual dipolar coupling `D` (both in Hz) from the peak
position difference between paired component spectra. `isotropic`/`aligned` are the plane
indices of the two components in each condition; `couplingdim` is the dimension the splitting
is measured in; `scale` is the fraction of the coupling the measured separation represents
(1 for IPAP, 0.5 for HSQC/TROSY); `gammasign` flips the sign for a negative-γ nucleus (¹⁵N).
"""
struct RDCModel <: FittingModel
    isotropic::Tuple{Int,Int}
    aligned::Tuple{Int,Int}
    couplingdim::Type
    scale::Float64
    gammasign::Float64
end

"""
    rdc2d(; isotropic, aligned, coupling=nothing, scale=1)

Interactive measurement of one-bond couplings and residual dipolar couplings from paired
2D spectra. Supply the two doublet-component spectra (already combined, e.g. IPAP α/β, or
HSQC/TROSY) for each condition:

    rdc2d(isotropic = ["iso_a/pdata/1", "iso_b/pdata/1"],
          aligned   = ["aln_a/pdata/1", "aln_b/pdata/1"])

The four spectra become the planes of a moving-peak experiment. Track each residue's peak
across the planes (T) or add and adjust by hand (A); the per-residue postfit then reports

    J   = sep(isotropic) / scale
    J+D = sep(aligned)   / scale
    D   = (J+D) − J

where `sep` is the position difference between the two components in the coupling dimension
(converted to Hz). `scale` is 1 for IPAP and 0.5 for HSQC/TROSY (the separation is then half
the coupling). `coupling` selects the dimension (`:F1`/`:F2`); it defaults to the
heteronuclear dimension, and the sign is flipped automatically for ¹⁵N (so J ≈ −93 Hz). List
the two components in the same order for both conditions; if J comes out with the wrong sign,
swap the pair.
"""
function rdc2d(; isotropic, aligned, coupling=nothing, scale=1.0)
    length(isotropic) == 2 || error("`isotropic` must be two component spectra, e.g. [A, B]")
    length(aligned) == 2 || error("`aligned` must be two component spectra, e.g. [A, B]")

    files = [isotropic[1], isotropic[2], aligned[1], aligned[2]]
    specdata = preparespecdata(files, IntensityExperiment)
    length(specdata.z) == 4 ||
        error("expected 4 single-plane spectra (got $(length(specdata.z))); each input must be a 2D spectrum")

    couplingdim = _resolve_coupling_dim(specdata, coupling)
    model = RDCModel((1, 2), (3, 4), couplingdim, Float64(scale),
                     _gamma_sign(specdata, couplingdim))

    peaks = Observable(Vector{Peak}())
    expt = MovingExperiment(specdata, peaks, model, nothing, CrossSectionVisualisation())
    return gui!(expt)
end

# Coupling dimension: :F1/:F2 (or :x/:y); default is the non-¹H (heteronuclear) dimension.
function _resolve_coupling_dim(specdata, coupling)
    (coupling === :F1 || coupling === :x) && return F1Dim
    (coupling === :F2 || coupling === :y) && return F2Dim
    isnothing(coupling) || error("`coupling` must be :F1/:F2 (or :x/:y), got $coupling")
    l1 = uppercase(string(label(specdata.nmrdata[1], F1Dim)))
    return occursin("H", l1) ? F2Dim : F1Dim
end

# Sign factor for the coupling: negative for ¹⁵N (negative gyromagnetic ratio), else positive.
function _gamma_sign(specdata, dim)
    return occursin("N", uppercase(string(label(specdata.nmrdata[1], dim)))) ? -1.0 : 1.0
end

function setup_post_parameters!(peak::Peak, ::RDCModel)
    peak.postparameters[:J] = Parameter("J", 0.0)
    peak.postparameters[:D] = Parameter("D", 0.0)
    return
end

function postfitmoving!(peak::Peak, expt::MovingExperiment, model::RDCModel)
    sym = model.couplingdim === F1Dim ? :x : :y
    axis = dims(expt.specdata.nmrdata[1], model.couplingdim)
    pos = peak.parameters[sym].value[]
    σ = peak.parameters[sym].uncertainty[]
    i1, i2 = model.isotropic
    a1, a2 = model.aligned

    # signed separations in Hz, and the ppm→Hz slope for converting position uncertainties
    sepiso = hz(pos[i2], axis) - hz(pos[i1], axis)
    sepaln = hz(pos[a2], axis) - hz(pos[a1], axis)
    slope = abs(hz(pos[i1] + 1.0, axis) - hz(pos[i1], axis))

    s = model.gammasign / model.scale
    J = s * sepiso
    D = s * (sepaln - sepiso)

    σiso = slope * sqrt(σ[i1]^2 + σ[i2]^2)
    σaln = slope * sqrt(σ[a1]^2 + σ[a2]^2)
    f = abs(s)
    peak.postparameters[:J].value[] .= J
    peak.postparameters[:J].uncertainty[] .= f * σiso
    peak.postparameters[:D].value[] .= D
    peak.postparameters[:D].uncertainty[] .= f * sqrt(σiso^2 + σaln^2)
    peak.postfitted[] = true
    return
end

primaryparam(expt::MovingExperiment) = expt.model isa RDCModel ? :D : :amp

function addpeakhint(::MovingPeakExperiment)
    return "Press (A) to add a peak, or (T) to add and track it across planes, under the cursor"
end

function peakinfotext(expt::MovingExperiment, idx)
    idx == 0 && return "No peak selected"
    peak = expt.peaks[][idx]
    if !peak.postfitted[]
        return "Peak: $(peak.label[])\nNot fitted"
    end

    if expt.model isa RDCModel
        J = peak.postparameters[:J]
        D = peak.postparameters[:D]
        return join(["Peak: $(peak.label[])",
                     "",
                     "J: $(round(J.value[][1]; digits=2)) ± $(round(J.uncertainty[][1]; digits=2)) Hz",
                     "D: $(round(D.value[][1]; digits=2)) ± $(round(D.uncertainty[][1]; digits=2)) Hz"],
                    "\n")
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
    # Above the contours/context, but below the drag handle (z=10) so the handle stays pickable.
    translate!(g[:plttrajectories], 0, 0, 2)
    translate!(g[:plttrajectorypts], 0, 0, 2)

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
    expt.model isa RDCModel &&
        return _rdc_summaryplot(expt; title=title, size=size,
                                include_unassigned=include_unassigned)

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

"""RDC summary: residual dipolar coupling D (Hz) against residue number."""
function _rdc_summaryplot(expt::MovingExperiment; title="", size=nothing,
                          include_unassigned=false)
    peaks = sortedpeaks(expt)
    if !include_unassigned
        assigned = filter(p -> extract_residue_number(p.label[]) > 0, peaks)
        isempty(assigned) || (peaks = assigned)
    end
    peaks = filter(p -> p.postfitted[] && haskey(p.postparameters, :D), peaks)

    figkw = isnothing(size) ? NamedTuple() : (; size=size)
    fig = Figure(; figkw...)
    ax = Axis(fig[1, 1]; xlabel="Residue number", ylabel="D / Hz", title=title,
              xgridvisible=false, ygridvisible=false)

    res = Float64[extract_residue_number(p.label[]) for p in peaks]
    D = Float64[p.postparameters[:D].value[][1] for p in peaks]
    Derr = Float64[p.postparameters[:D].uncertainty[][1] for p in peaks]
    errorbars!(ax, res, D, Derr; whiskerwidth=6, color=:black)
    scatter!(ax, res, D; color=:steelblue)
    hlines!(ax, [0]; linewidth=0)
    return fig
end
