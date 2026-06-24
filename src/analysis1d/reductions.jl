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
    noiselevel(trace, noise) -> Float64

Per-point RMS noise estimated as the standard deviation of intensities in the noise
region.
"""
noiselevel(t::Trace, noise::Region) = std(@view t.y[roi_indices(t, noise)])

"""
    reduce_region(::Integrate, region, dataset) -> NamedTuple

Integrate `region` over every plane of `dataset`, returning a `NamedTuple` of named
quantity series. For `Integrate` there is a single series `I`, a
`Vector{Measurement}` whose uncertainties propagate the per-point noise over the
number of summed points (σ_I = σ_point · √N_roi). Returning a `NamedTuple` keeps the
contract general: a future NMF reduction can return several named component series.
"""
function reduce_region(::Integrate, region::Region, dataset::Dataset1D)
    planes = dataset.planes
    I = map(planes.traces) do t
        val = integrate(t, region)
        n = length(roi_indices(t, region))
        σ = noiselevel(t, dataset.noise) * sqrt(n)
        return val ± σ
    end
    return (; I)
end

"""
    integrals(region, dataset) -> Vector{Measurement}

Convenience accessor for the primary (`I`) series of an `Integrate` reduction.
"""
integrals(region::Region, dataset::Dataset1D) = reduce_region(Integrate(), region, dataset).I
