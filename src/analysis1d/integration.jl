"""
    roiindices(region, trace) -> Vector{Int}

Indices of `trace.δ` falling within `region`. For a zero-width region (or one that
contains no grid points) the single nearest point is returned, giving a peak height.
"""
function roiindices(r::Region, t::Trace)
    idx = findall(δ -> r.lo ≤ δ ≤ r.hi, t.δ)
    isempty(idx) || return idx
    mid = (r.lo + r.hi) / 2
    return [argmin(abs.(t.δ .- mid))]
end

"""
    integrate(region, trace) -> Float64

Summed intensity of `region` in one plane (a height when the region has zero width — the
nearest point is taken). The summed intensity is the convention used by the legacy 1D
routines and by GUI2D's intensity analysis.
"""
integrate(r::Region, t::Trace) = sum(@view t.y[roiindices(r, t)])

"""
    integrate(region, dataset) -> Vector{Measurement{Float64}}

Integrate `region` over every plane of `dataset`, with an uncertainty attached to each.

The uncertainty is measured rather than taken from `spec[:noise]`: the noise on an
*integrated* region depends on its width and on the noise where it actually sits, not on
the whole-spectrum RMS level. So a noise region of the same width, centred at
`dataset.noisecenter`, is integrated over every plane and the standard deviation of those
integrals taken. With uniform noise this approaches `sqrt(n points)`, the trace
intensities already being noise-normalised, but it tracks a locally noisier baseline.
"""
function integrate(region::Region, ds::Dataset1D)
    w = width(region)
    w == 0 && (w = 0.05)   # a height still needs a nominal window for the noise estimate
    noiseregion = Region("noise", ds.noisecenter - w / 2, ds.noisecenter + w / 2)

    raw = [integrate(region, t) for t in ds.planes.traces]
    noiseintegrals = [integrate(noiseregion, t) for t in ds.planes.traces]
    σ = length(noiseintegrals) > 1 ? std(noiseintegrals) : abs(noiseintegrals[1])
    (σ == 0 || isnan(σ)) && (σ = 1.0)

    return [v ± σ for v in raw]
end
