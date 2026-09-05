"""
    Reduction

Abstract supertype for the rule mapping a region across all planes to a series of
quantities (one or more per plane). v1 provides only [`Integrate`](@ref). Future
reductions (e.g. NMF for kinetics) plug in here by implementing
[`reduceregion`](@ref).
"""
abstract type Reduction end

"""
    Integrate()

Reduce a region to its integrated intensity per plane (a height when the region has
zero width — the nearest point is taken). The summed intensity is the convention used
by the existing 1D routines and by GUI2D's intensity analysis.
"""
struct Integrate <: Reduction end

"""
    roiindices(trace, region) -> Vector{Int}

Indices of `trace.δ` falling within `region`. For a zero-width region (or one that
contains no grid points) the single nearest point is returned, giving a peak height.
"""
function roiindices(t::Trace, r::Region)
    idx = findall(δ -> r.lo ≤ δ ≤ r.hi, t.δ)
    isempty(idx) || return idx
    mid = (r.lo + r.hi) / 2
    return [argmin(abs.(t.δ .- mid))]
end

"""
    integrate(trace, region) -> Float64

Summed intensity of `trace` over `region`.
"""
integrate(t::Trace, r::Region) = sum(@view t.y[roiindices(t, r)])

"""
    reduceregion(::Integrate, region, dataset) -> NamedTuple

Integrate `region` over every plane of `dataset`, returning a `NamedTuple` of named
quantity series. For `Integrate` there is a single series `I`, a `Vector{Measurement}`.

This is a different noise question from the point-wise spectrum normalisation done in
`tracesfromspec` (dividing by `spec[:noise]`, the whole-spectrum RMS level): the
uncertainty on an *integrated* region depends on its width and on the actual, possibly
non-white noise there, not on a single global scalar. So, following the legacy 1D
routines and `Exchange1D`, it's measured directly: integrate a noise region of the same
width as `region` (centred at `dataset.noisecenter`, which the GUI lets the user
reposition) over every plane, and take the standard deviation of those integrals across
planes. Since the trace intensities are already noise-normalised, this comes out close
to `sqrt(n points)` when the noise is uniform, but tracks the real, locally-measured
noise otherwise. Returning a `NamedTuple` keeps the contract general: a future NMF
reduction can return several named component series.
"""
function reduceregion(::Integrate, region::Region, dataset::Dataset1D)
    planes = dataset.planes
    w = width(region)
    w == 0 && (w = 0.05)   # a height still needs a nominal window for the noise estimate
    noiseregion = Region("noise", dataset.noisecenter - w / 2, dataset.noisecenter + w / 2)

    raw = [integrate(t, region) for t in planes.traces]
    noiseintegrals = [integrate(t, noiseregion) for t in planes.traces]
    σ = length(noiseintegrals) > 1 ? std(noiseintegrals) : abs(noiseintegrals[1])
    (σ == 0 || isnan(σ)) && (σ = 1.0)

    I = [v ± σ for v in raw]
    return (; I)
end

"""
    integrals(region, dataset) -> Vector{Measurement}

Convenience accessor for the primary (`I`) series of an `Integrate` reduction.
"""
integrals(region::Region, dataset::Dataset1D) = reduceregion(Integrate(), region, dataset).I
