"""
    Reduction

Abstract supertype for the rule mapping a region across all planes to a series of
quantities (one or more per plane). v1 provides only [`Integrate`](@ref). Future
reductions (e.g. NMF for kinetics) plug in here by implementing
[`reduce_region`](@ref).
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
    roi_indices(trace, region) -> Vector{Int}

Indices of `trace.δ` falling within `region`. For a zero-width region (or one that
contains no grid points) the single nearest point is returned, giving a peak height.
"""
function roi_indices(t::Trace, r::Region)
    idx = findall(δ -> r.lo ≤ δ ≤ r.hi, t.δ)
    isempty(idx) || return idx
    mid = (r.lo + r.hi) / 2
    return [argmin(abs.(t.δ .- mid))]
end

"""
    integrate(trace, region) -> Float64

Summed intensity of `trace` over `region`.
"""
integrate(t::Trace, r::Region) = sum(@view t.y[roi_indices(t, r)])

"""
    reduce_region(::Integrate, region, dataset) -> NamedTuple

Integrate `region` over every plane of `dataset`, returning a `NamedTuple` of named
quantity series. For `Integrate` there is a single series `I`, a `Vector{Measurement}`.

Following the legacy 1D routines and `Exchange1D`, the noise estimate is the standard
deviation, taken *across planes*, of the integral over a noise region of the **same
width** as `region` (centred at `dataset.noisecenter`) — a region of matching width
integrates the same amount of noise as the signal region, so this directly estimates the
per-plane integration noise without needing a separate width control. Intensities are
then scaled down by that noise level, so values read in signal-to-noise units instead of
raw (often very large) integrated intensities; every point in the series then carries
uncertainty exactly 1 by construction. Returning a `NamedTuple` keeps the contract
general: a future NMF reduction can return several named component series.
"""
function reduce_region(::Integrate, region::Region, dataset::Dataset1D)
    planes = dataset.planes
    w = width(region)
    w == 0 && (w = 0.05)   # a height still needs a nominal window for the noise estimate
    noiseregion = Region("noise", dataset.noisecenter - w / 2, dataset.noisecenter + w / 2)

    raw = [integrate(t, region) for t in planes.traces]
    noiseintegrals = [integrate(t, noiseregion) for t in planes.traces]
    σ = length(noiseintegrals) > 1 ? std(noiseintegrals) : abs(noiseintegrals[1])
    (σ == 0 || isnan(σ)) && (σ = 1.0)

    I = [(v / σ) ± 1.0 for v in raw]
    return (; I)
end

"""
    integrals(region, dataset) -> Vector{Measurement}

Convenience accessor for the primary (`I`) series of an `Integrate` reduction.
"""
integrals(region::Region, dataset::Dataset1D) = reduce_region(Integrate(), region, dataset).I
