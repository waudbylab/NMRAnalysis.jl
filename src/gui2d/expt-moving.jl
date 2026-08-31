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

# Moving-peak spectra store per-plane S/N (data/noise) so contour levels are consistent across
# heterogeneous planes (see preparespecdata). For cross-sections we instead want intensities
# that are comparable across planes, so undo the per-plane noise normalisation and rescale by
# scale(spec) (number of scans × receiver gain): data/noise × noise/scale = data/scale.
function crosssection_scale(expt::MovingExperiment, i)
    spec = expt.specdata.nmrdata[i]
    return spec[:noise] / scale(spec)
end

# Whether the (T) add-and-track command applies. Tracking the intensity maximum makes sense
# for a single peak walking across planes (titrations) but not for paired component spectra,
# so it is disabled for RDC experiments.
cantrack(::Experiment) = false
cantrack(expt::MovingExperiment) = !(expt.model isa RDCModel)

# Whether adding a peak should walk plane-by-plane (the default) or place immediately using the
# learned displacement pattern. Generic moving peaks always walk; RDC only walks until a peak
# has been fitted, after which new peaks copy the established offsets (addpeak! seeds from
# average_displacements) so the user need not re-mark every plane.
needsguidedadd(expt::MovingExperiment) =
    !(expt.model isa RDCModel) || !any(p -> p.postfitted[], expt.peaks[])

# Moving-peak spectra are normalised per plane by each plane's OWN noise level, so contour
# levels and signal-to-noise are consistent even when the planes are different experiment
# types (e.g. rdc2d's isotropic vs aligned spectra, or HSQC vs TROSY). This differs from the
# intensity experiments (relaxation, CEST, ...), which normalise every plane by the first
# plane's noise to keep amplitudes directly comparable across the series.
function preparespecdata(inputfilenames, ::Type{MovingExperiment})
    @debug "Preparing spec data for moving-peak experiment: $inputfilenames"
    spec, x, y, z, σ, zlabels = if inputfilenames isa String
        spec, x, y, z, σ = loadspecdata(inputfilenames, IntensityExperiment)
        (SingleElementVector(spec),
         SingleElementVector(x),
         SingleElementVector(y),
         z ./ σ,
         SingleElementVector(1),
         SingleElementVector(choptitle(label(spec))))
    elseif inputfilenames isa Vector{String}
        tmp = loadspecdata.(inputfilenames, IntensityExperiment)
        spec = []
        x = []
        y = []
        z = []
        σ = []
        zlabels = []
        for t in tmp
            n = length(t[4])
            if n == 1
                push!(spec, t[1])
                push!(x, t[2])
                push!(y, t[3])
                push!(z, t[4][1] ./ t[5])  # normalise by THIS plane's noise
                push!(σ, t[5])
                push!(zlabels, choptitle(label(t[1])))
            else
                append!(spec, fill(t[1], n))
                append!(x, fill(t[2], n))
                append!(y, fill(t[3], n))
                append!(z, [zi ./ t[5] for zi in t[4]])
                append!(σ, fill(t[5], n))
                append!(zlabels, fill(choptitle(label(t[1])), n))
            end
        end
        map(MaybeVector, (spec, x, y, z, σ, zlabels))
    end
    return SpecData(spec, x, y, z, σ, zlabels)
end

"""
    peaktrack2d(inputfilenames)

Start an interactive GUI for analysing a series of 2D spectra in which **peak positions
change** from plane to plane (e.g. a titration, or a coupling-constant / RDC measurement).
Each peak is fitted to a 2D Lorentzian lineshape independently in every plane, so its
position and linewidth are free to move; no physical model is applied to the trajectory.

Use this to measure how peak positions and linewidths evolve across a series for downstream
analysis, or as the basis for the position-based physical models (titration, coupling).

# Arguments
- `inputfilenames`: A single path string (pseudo-3D dataset) or a vector of path strings
  (one file per plane) pointing to processed Bruker data directories.

# Example
```julia
peaktrack2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"])
```

See also [`titration2d`](@ref) for fitting binding isotherms to a titration series.
"""
function peaktrack2d(inputfilenames)
    specdata = preparespecdata(inputfilenames, MovingExperiment)
    peaks = Observable(Vector{Peak}())

    expt = MovingExperiment(specdata, peaks, NoFitting(), nothing)

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

    # Seed the other planes from the average per-plane displacement of already-fitted peaks
    # (relative to the plane being clicked in), so a new peak inherits the common motion pattern
    # - e.g. the IPAP/RDC zig-zag - instead of starting flat at the click. Falls back to the flat
    # seed for the first peak. (T) overrides this with maximum-tracking.
    disp = average_displacements(expt)
    if !isnothing(disp)
        dx, dy = disp
        for i in 1:n
            setpeakposition!(expt, newpeak, i, x0 + dx[i], y0 + dy[i])
        end
    end

    # Post-parameters based on model type (none for NoFitting)
    setup_post_parameters!(newpeak, expt.model)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""
    average_displacements(expt) -> (dx, dy) | nothing

Mean per-plane position displacement of the already-fitted peaks, relative to the current
plane (so `dx`/`dy` are zero at that plane). Used to seed a new peak's planes with the common
motion pattern. Returns `nothing` if no peak has been fitted yet.
"""
function average_displacements(expt::MovingExperiment)
    fitted = [p for p in expt.peaks[] if p.postfitted[]]
    isempty(fitted) && return nothing
    n = nslices(expt)
    anchor = expt.state[][:current_slice][]
    dx = zeros(n)
    dy = zeros(n)
    for p in fitted
        xv = p.parameters[:x].value[]
        yv = p.parameters[:y].value[]
        for i in 1:n
            dx[i] += xv[i] - xv[anchor]
            dy[i] += yv[i] - yv[anchor]
        end
    end
    return dx ./ length(fitted), dy ./ length(fitted)
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

# --- guided add (A) ----------------------------------------------------------
# Adding a moving peak walks through the planes so the user marks its position in each one,
# starting at the current plane and stepping forwards (wrapping back to the start), then
# returning to the original plane. While adding (mode :adding): (a) marks this plane and
# advances, (space) fills the remaining planes with the current position and finishes, (esc)
# cancels. State is held in the GUI state dict under :add_*.

# Publish the marked-so-far positions for the in-progress add overlay (see add_moving_overlays!).
function show_add_marks!(state, positions)
    haskey(state, :add_marks) || return
    return state[:add_marks][] = [p for p in positions if !isnan(p[1])]
end

"""Begin a guided add at the current plane, marking it immediately at `pos`."""
function begin_add!(expt::MovingPeakExperiment, state, pos)
    n = nslices(expt)
    s = state[:current_slice][]
    state[:add_order] = Observable([mod1(s + k, n) for k in 0:(n - 1)])  # current, then forwards
    state[:add_index] = Observable(1)
    state[:add_positions] = Observable(fill(Point2f(NaN32, NaN32), n))
    state[:add_anchor] = Observable(s)
    state[:mode][] = :adding
    return mark_add!(expt, state, pos)
end

"""Record `pos` for the plane currently being marked and advance (or finish)."""
function mark_add!(expt::MovingPeakExperiment, state, pos)
    n = nslices(expt)
    order = state[:add_order][]
    idx = state[:add_index][]
    positions = copy(state[:add_positions][])
    positions[order[idx]] = Point2f(pos)
    state[:add_positions][] = positions
    show_add_marks!(state, positions)
    if idx >= n
        return finish_add!(expt, state)
    end
    state[:add_index][] = idx + 1
    return set_close_to!(state[:gui][][:sliderslice], order[idx + 1])
end

"""Fill every not-yet-marked plane with `pos`, then finish."""
function fill_add!(expt::MovingPeakExperiment, state, pos)
    positions = copy(state[:add_positions][])
    for i in eachindex(positions)
        isnan(positions[i][1]) && (positions[i] = Point2f(pos))
    end
    state[:add_positions][] = positions
    show_add_marks!(state, positions)
    return finish_add!(expt, state)
end

"""Create the peak from the marked positions, return to the original plane, leave :adding."""
function finish_add!(expt::MovingPeakExperiment, state)
    positions = state[:add_positions][]
    anchor = state[:add_anchor][]
    fallback = positions[anchor]
    addpeak!(expt, Point2f(fallback))
    peak = expt.peaks[][end]
    for i in eachindex(positions)
        p = isnan(positions[i][1]) ? fallback : positions[i]
        setpeakposition!(expt, peak, i, p[1], p[2])
    end
    peak.touched[] = true
    notify(expt.peaks)
    show_add_marks!(state, Point2f[])
    state[:current_peak_idx][] = length(expt.peaks[])
    state[:mode][] = :normal
    return set_close_to!(state[:gui][][:sliderslice], anchor)
end

"""Abort the guided add and return to the original plane."""
function cancel_add!(expt::MovingPeakExperiment, state)
    anchor = state[:add_anchor][]
    show_add_marks!(state, Point2f[])
    state[:mode][] = :normal
    return set_close_to!(state[:gui][][:sliderslice], anchor)
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

The four spectra become the planes of a peak-tracking experiment. Track each residue's peak
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
    specdata = preparespecdata(files, MovingExperiment)
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

# --- titration / binding isotherm -------------------------------------------

"""
    TitrationModel

Fit a single-site (1:1) binding isotherm to a titration series, deriving a **global**
dissociation constant `Kd` (shared by all residues) together with per-residue, per-dimension
free and bound chemical shifts (`δfree`, `δbound`).

The bound fraction in each plane is computed from the ligand concentration (`expt.x`) and, if
available, the protein concentration (`protein`, one per plane). With protein concentrations
the exact 1:1 quadratic is used; with ligand concentrations only, the hyperbolic
`f = [L]/(Kd + [L])` approximation (ligand ≈ free) is used. The two are consistent in the
limit `[P] → 0`.

# Fields
- `protein`: Per-plane total protein concentration, or `nothing` when only ligand
  concentrations were supplied.
- `weights`: `(wx, wy)` dimension weights for the combined CSP (Williamson 2013); the default
  `(1.0, 0.14)` assumes the F1 (x) axis is ¹H and the F2 (y) axis is ¹⁵N.
"""
struct TitrationModel <: FittingModel
    protein::Union{Vector{Float64},Nothing}
    weights::Tuple{Float64,Float64}
end

"""
    titration2d(inputfilenames; L0, P0=nothing, weights=(1.0, 0.14))

Interactive analysis of a 2D titration series, fitting a global binding isotherm to the
chemical-shift perturbations. Built on [`peaktrack2d`](@ref): track each residue's peak across
the planes (`T`) or add and adjust by hand (`A`), and the fit recovers a single `Kd` shared by
all residues plus per-residue free/bound shifts in each dimension.

# Arguments
- `inputfilenames`: A single path string (pseudo-3D dataset) or a vector of path strings (one
  file per plane) pointing to processed Bruker data directories.
- `L0`: Total **ligand** concentration in each plane (one per plane).
- `P0`: Total **protein** concentration in each plane (one per plane), or `nothing`. When given,
  the exact 1:1 binding equation is used, accounting for protein concentration and dilution;
  otherwise the hyperbolic (ligand ≈ free) approximation is used. `Kd` is reported in the same
  concentration units as `L0`/`P0`.
- `weights`: `(wx, wy)` weighting of the two dimensions for the combined CSP `|Δδ|`; the default
  assumes ¹H (x) / ¹⁵N (y).

# Results
- The per-residue panel shows ΔδX and ΔδN against ligand concentration (referenced to the
  fitted `δfree`), with the fitted binding curve overlaid.
- The global `Kd` is reported in the results panel.
- The summary plot shows ΔδX, ΔδN and the combined `|Δδ|` (saturation CSP = `δbound − δfree`)
  against residue number.

# Example
```julia
titration2d(files; L0=ligand_concs)                  # ligand concentrations only
titration2d(files; L0=ligand_concs, P0=protein_concs) # exact 1:1, accounts for dilution
```
"""
function titration2d(inputfilenames; L0, P0=nothing, weights=(1.0, 0.14))
    isempty(L0) && error("`L0` (ligand concentrations) must not be empty")

    ligand = Float64.(collect(L0))
    protein = isnothing(P0) ? nothing : Float64.(collect(P0))

    specdata = preparespecdata(inputfilenames, MovingExperiment)
    length(specdata.z) == length(ligand) ||
        error("got $(length(ligand)) ligand concentrations (L0) for $(length(specdata.z)) planes")
    isnothing(protein) || length(protein) == length(ligand) ||
        error("got $(length(protein)) protein concentrations (P0) for $(length(ligand)) ligand concentrations (L0)")

    model = TitrationModel(protein, (Float64(weights[1]), Float64(weights[2])))
    peaks = Observable(Vector{Peak}())
    expt = MovingExperiment(specdata, peaks, model, ligand, TitrationVisualisation())
    return gui!(expt)
end

# Fraction of protein bound for total protein `Pt` (or `nothing`) and total ligand `L`. With a
# protein concentration the exact 1:1 quadratic is used; otherwise the hyperbolic (ligand ≈
# free) approximation. Both reduce to the same form as [P] → 0.
function _fraction(Kd, Pt, L)
    Kd = max(Kd, 1e-12)
    (isnothing(Pt) || Pt ≤ 0) && return L / (Kd + L)
    s = Pt + L + Kd
    pl = (s - sqrt(max(s^2 - 4 * Pt * L, 0.0))) / 2
    return pl / Pt
end

# Bound fraction in plane `i` of a titration experiment.
boundfraction(model::TitrationModel, Lt, i, Kd) =
    _fraction(Kd, isnothing(model.protein) ? nothing : model.protein[i], Lt[i])

# Weighted linear regression of `y` on `f`: returns (a, b) with model y ≈ a + b·f. `w` are
# per-point weights. Used both inside the Kd search (variable projection) and to read off the
# per-residue free shift (a = δfree) and span (b = δbound − δfree).
function _wlinreg(f, y, w)
    sw = sum(w)
    sf = sum(w .* f)
    sy = sum(w .* y)
    sff = sum(w .* f .* f)
    sfy = sum(w .* f .* y)
    det = sw * sff - sf^2
    abs(det) < 1e-30 && return (sy / sw, 0.0)
    b = (sw * sfy - sf * sy) / det
    a = (sy - b * sf) / sw
    return (a, b)
end

# Per-point weights for one residue/dimension series: inverse position uncertainty where the
# per-plane fit produced a finite one, falling back to the series' median (or unity) so a few
# bad σ don't dominate or zero out the series.
function _seriesweights(σ)
    good = [s for s in σ if isfinite(s) && s > 0]
    fallback = isempty(good) ? 1.0 : median(good)
    return [isfinite(s) && s > 0 ? 1.0 / s : 1.0 / fallback for s in σ]
end

setup_post_parameters!(peak::Peak, ::TitrationModel) = setup_titration_postparameters!(peak)

function setup_titration_postparameters!(peak::Peak)
    peak.postparameters[:CSP] = Parameter("CSP", 0.0)      # combined |Δδ| (weighted)
    peak.postparameters[:dX] = Parameter("dX", 0.0)        # ΔδX = δbound − δfree (x dim)
    peak.postparameters[:dY] = Parameter("dY", 0.0)        # ΔδY = δbound − δfree (y dim)
    peak.postparameters[:Xfree] = Parameter("Xfree", 0.0)
    peak.postparameters[:Xbound] = Parameter("Xbound", 0.0)
    peak.postparameters[:Yfree] = Parameter("Yfree", 0.0)
    peak.postparameters[:Ybound] = Parameter("Ybound", 0.0)
    peak.postparameters[:Kd] = Parameter("Kd", 0.0)
    return
end

# Per-peak postfit is a no-op beyond flagging the peak fitted: the binding isotherm is fitted
# globally across all residues (see postfitglobal!).
function postfitmoving!(peak::Peak, ::MovingExperiment, ::TitrationModel)
    peak.postfitted[] = true
    return
end

# Global binding-isotherm fit. The only nonlinear parameter is the shared Kd: for any Kd the
# per-residue/dimension δfree and δbound are obtained by (weighted) linear regression of the
# fitted positions on the bound fraction (variable projection). So we optimise Kd over a 1-D
# residual built from all residues and dimensions, then read off the linear parameters.
postfitglobal!(expt::MovingExperiment) = postfitglobal!(expt, expt.model)
postfitglobal!(::MovingExperiment, ::FittingModel) = nothing

function postfitglobal!(expt::MovingExperiment, model::TitrationModel)
    peaks = [p for p in expt.peaks[] if p.postfitted[]]
    isempty(peaks) && return
    Lt = expt.x
    n = nslices(expt)

    # one (measured, weights) series per residue per dimension
    series = NamedTuple{(:y, :w),Tuple{Vector{Float64},Vector{Float64}}}[]
    for peak in peaks, sym in (:x, :y)
        y = Float64.(peak.parameters[sym].value[])
        w = _seriesweights(peak.parameters[sym].uncertainty[])
        push!(series, (y=y, w=w))
    end

    fractions(Kd) = [boundfraction(model, Lt, i, Kd) for i in 1:n]

    # residual vector across all series at a given Kd (linear params projected out)
    function resid(p)
        Kd = exp(p[1])
        f = fractions(Kd)
        r = Float64[]
        for s in series
            a, b = _wlinreg(f, s.y, s.w)
            append!(r, sqrt.(s.w) .* (s.y .- (a .+ b .* f)))
        end
        return r
    end

    Kd0 = _initial_kd(Lt)
    local Kd, Kderr
    try
        sol = LsqFit.lmfit(resid, [log(Kd0)], Float64[]; autodiff=:finite, maxIter=200)
        logKd = coef(sol)[1]
        Kd = exp(logKd)
        Kderr = Kd * stderror(sol)[1]   # delta method: σ(Kd) = Kd·σ(logKd)
    catch e
        e isa FitCancelled && rethrow()
        @debug "Titration global fit failed" exception = e
        Kd = Kd0
        Kderr = Inf
    end

    # read off per-residue/dimension free & bound shifts at the fitted Kd
    f = fractions(Kd)
    wx, wy = model.weights
    si = 1
    for peak in peaks
        ax, bx = _wlinreg(f, series[si].y, series[si].w)
        ay, by = _wlinreg(f, series[si + 1].y, series[si + 1].w)
        si += 2

        peak.postparameters[:Xfree].value[] .= ax
        peak.postparameters[:Xbound].value[] .= ax + bx
        peak.postparameters[:Yfree].value[] .= ay
        peak.postparameters[:Ybound].value[] .= ay + by
        peak.postparameters[:dX].value[] .= bx
        peak.postparameters[:dY].value[] .= by
        peak.postparameters[:CSP].value[] .= sqrt((wx * bx)^2 + (wy * by)^2)
        peak.postparameters[:Kd].value[] .= Kd
        peak.postparameters[:Kd].uncertainty[] .= Kderr
    end
    return
end

# Initial Kd guess: the median non-zero ligand concentration (roughly the half-saturation
# point for a hyperbolic curve), falling back to 1.0.
function _initial_kd(Lt)
    nz = [l for l in Lt if l > 0]
    return isempty(nz) ? 1.0 : median(nz)
end

primaryparam(expt::MovingExperiment) =
    expt.model isa RDCModel ? :D : expt.model isa TitrationModel ? :CSP : :amp

function addpeakhint(expt::MovingPeakExperiment)
    s = "Press (A) to add a peak, marking its position in each plane"
    cantrack(expt) && (s *= ", or (T) to add and auto-track across planes")
    return s
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

    if expt.model isa TitrationModel
        Kd = peak.postparameters[:Kd]
        dX = peak.postparameters[:dX].value[][1]
        dY = peak.postparameters[:dY].value[][1]
        csp = peak.postparameters[:CSP].value[][1]
        return join(["Peak: $(peak.label[])",
                     "",
                     "Kd: $(round(Kd.value[][1]; digits=3)) ± $(round(Kd.uncertainty[][1]; digits=3)) (global)",
                     "",
                     "ΔδX: $(round(dX; digits=4)) ppm",
                     "ΔδY: $(round(dY; digits=4)) ppm",
                     "|Δδ|: $(round(csp; digits=4)) ppm"],
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
    analysis = expt.model isa TitrationModel ? "Titration" :
               expt.model isa RDCModel ? "RDC / coupling" : "Peak tracking"
    info = ["Analysis type: $analysis",
            "Model: $(typeof(expt.model))",
            "Filename: $(expt.specdata.nmrdata[1][:filename])",
            "Number of peaks: $(length(expt.peaks[]))",
            "Number of planes: $(nslices(expt))"]
    if expt.model isa TitrationModel
        fitted = [p for p in expt.peaks[] if p.postfitted[]]
        if !isempty(fitted)
            kd = fitted[1].postparameters[:Kd]
            push!(info,
                  "Global Kd: $(round(kd.value[][1]; digits=4)) ± $(round(kd.uncertainty[][1]; digits=4))")
        end
    end
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

    # In-progress guided-add marks: orange crosses at the positions marked so far.
    state[:add_marks] = Observable(Point2f[])
    g[:pltaddmarks] = scatter!(ax, state[:add_marks]; marker=:xcross, markersize=16,
                               color=:darkorange)
    translate!(g[:pltaddmarks], 0, 0, 11)

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
    expt.model isa TitrationModel &&
        return _titration_summaryplot(expt; title=title, size=size,
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

"""
Titration summary: three stacked panels of the fitted saturation perturbation against residue
number — ΔδX (= δbound − δfree, x dim), ΔδN (y dim, assuming ¹⁵N) and the weighted combination
|Δδ|. Values are the extrapolated full-saturation CSPs from the global fit, so the experiment
must have been fitted.
"""
function _titration_summaryplot(expt::MovingExperiment; title="", size=nothing,
                                include_unassigned=false)
    peaks = sortedpeaks(expt)
    if !include_unassigned
        assigned = filter(p -> extract_residue_number(p.label[]) > 0, peaks)
        isempty(assigned) || (peaks = assigned)
    end
    peaks = filter(p -> p.postfitted[] && haskey(p.postparameters, :CSP), peaks)

    figkw = isnothing(size) ? NamedTuple() : (; size=size)
    fig = Figure(; figkw...)

    res = Float64[extract_residue_number(p.label[]) for p in peaks]
    panels = ((:dX, "ΔδX / ppm"), (:dY, "ΔδN / ppm"), (:CSP, "|Δδ| / ppm"))
    axes = Axis[]
    for (row, (sym, ylab)) in enumerate(panels)
        ax = Axis(fig[row, 1]; xlabel=(row == 3 ? "Residue number" : ""), ylabel=ylab,
                  title=(row == 1 ? title : ""), xgridvisible=false, ygridvisible=false)
        vals = Float64[p.postparameters[sym].value[][1] for p in peaks]
        scatterlines!(ax, res, vals; color=:steelblue)
        hlines!(ax, [0]; linewidth=0)
        push!(axes, ax)
    end
    length(axes) > 1 && linkxaxes!(axes...)
    return fig
end

# --- titration per-peak panel (residue highlight) ---------------------------
# The selected residue's chemical-shift perturbation against ligand concentration, one panel
# per dimension (ΔδX and ΔδY), with the fitted binding curve overlaid. Δδ is referenced to the
# fitted δfree once the global fit has run, and to plane 1 beforehand. Δδ is signed (not |Δδ|).

struct TitrationVisualisation <: VisualisationStrategy end

# Clamped piecewise-linear interpolation of `ys` (sampled at `xs`) at query point `q`.
function _lininterp(xs, ys, q)
    q ≤ xs[1] && return ys[1]
    q ≥ xs[end] && return ys[end]
    k = searchsortedfirst(xs, q)
    x0, x1 = xs[k - 1], xs[k]
    y0, y1 = ys[k - 1], ys[k]
    x1 == x0 && return y0
    return y0 + (y1 - y0) * (q - x0) / (x1 - x0)
end

# Smooth binding curve Δδ(L) = f(L)·(δbound − δfree) over a dense ligand grid. When protein
# concentrations are present they are interpolated against ligand so the quadratic can be
# evaluated along the curve.
function _titration_curve(model::TitrationModel, Lt, Kd, δfree, δbound; npts=100)
    lo = min(0.0, minimum(Lt))
    hi = maximum(Lt)
    hi ≤ lo && (hi = lo + 1.0)
    grid = range(lo, hi + 0.05 * (hi - lo), npts)
    Pg = if isnothing(model.protein)
        nothing
    else
        o = sortperm(Lt)
        [_lininterp(Lt[o], model.protein[o], g) for g in grid]
    end
    return [Point2f(g, _fraction(Kd, isnothing(Pg) ? nothing : Pg[k], g) * (δbound - δfree))
            for (k, g) in enumerate(grid)]
end

"""
    get_titration_data(peak, expt) -> (obsX, obsY, fitX, fitY)

Per-dimension Δδ-vs-concentration points for the highlighted residue: observed points
(`obsX`/`obsY`) and the fitted binding curve (`fitX`/`fitY`, empty until the global fit runs).
"""
function get_titration_data(peak, expt::MovingExperiment)
    isnothing(peak) && return (Point2f[], Point2f[], Point2f[], Point2f[])
    Lt = expt.x
    n = nslices(expt)
    xv = peak.parameters[:x].value[]
    yv = peak.parameters[:y].value[]

    fitted = peak.postfitted[] && haskey(peak.postparameters, :Kd)
    δXfree = fitted ? peak.postparameters[:Xfree].value[][1] : xv[1]
    δYfree = fitted ? peak.postparameters[:Yfree].value[][1] : yv[1]

    obsX = [Point2f(Lt[i], xv[i] - δXfree) for i in 1:n]
    obsY = [Point2f(Lt[i], yv[i] - δYfree) for i in 1:n]

    if fitted
        Kd = peak.postparameters[:Kd].value[][1]
        fitX = _titration_curve(expt.model, Lt, Kd, δXfree,
                                peak.postparameters[:Xbound].value[][1])
        fitY = _titration_curve(expt.model, Lt, Kd, δYfree,
                                peak.postparameters[:Ybound].value[][1])
    else
        fitX = Point2f[]
        fitY = Point2f[]
    end
    return (obsX, obsY, fitX, fitY)
end

function completestate!(state, expt::MovingExperiment, ::TitrationVisualisation)
    state[:peak_plot_data] = lift(peak -> get_titration_data(peak, expt), state[:current_peak])
    state[:peak_plot_obsX] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_obsY] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_fitX] = lift(d -> d[3], state[:peak_plot_data])
    return state[:peak_plot_fitY] = lift(d -> d[4], state[:peak_plot_data])
end

function makepeakplot!(gui, state, expt::MovingExperiment, ::TitrationVisualisation)
    gui[:axpeakplotX] = axX = Axis(gui[:panelpeakplot][1, 1];
                                   xlabel="[ligand]", ylabel="ΔδX / ppm")
    gui[:axpeakplotY] = axY = Axis(gui[:panelpeakplot][1, 2];
                                   xlabel="[ligand]", ylabel="ΔδY / ppm")
    hlines!(axX, [0]; linewidth=0)
    lines!(axX, state[:peak_plot_fitX]; color=:red)
    scatter!(axX, state[:peak_plot_obsX])
    hlines!(axY, [0]; linewidth=0)
    lines!(axY, state[:peak_plot_fitY]; color=:red)
    return scatter!(axY, state[:peak_plot_obsY])
end

function plot_peak!(panel, peak, expt::MovingExperiment, ::TitrationVisualisation)
    obsX, obsY, fitX, fitY = get_titration_data(peak, expt)
    axX = Axis(panel[1, 1]; xlabel="[ligand]", ylabel="ΔδX / ppm")
    axY = Axis(panel[1, 2]; xlabel="[ligand]", ylabel="ΔδY / ppm")
    hlines!(axX, [0]; linewidth=0)
    lines!(axX, fitX; color=:red)
    scatter!(axX, obsX)
    hlines!(axY, [0]; linewidth=0)
    lines!(axY, fitY; color=:red)
    return scatter!(axY, obsY)
end
