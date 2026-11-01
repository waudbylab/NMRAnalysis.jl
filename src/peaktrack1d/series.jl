# The spectra behind the stack, as plain arrays on one shared chemical shift axis.

"""
    SpectrumSeries(shifts, intensities, noise, sources, experiments; nucleus="1H", title="")
    SpectrumSeries(spectra)

The data behind the stack: a shared chemical shift axis (ppm), the real intensities with
one column per plane, the rms noise of each plane, the file each plane came from, and the
index of the experiment each plane belongs to. A pseudo-2D experiment contributes several
planes with one experiment index, which is what lets it be processed with one phase.

Built from a vector of processed spectra, 1D or pseudo-2D, whose shift axes must agree.
Intensities and noise are divided by NMRTools' `scale(spec)` (number of scans × receiver
gain, and the sample concentration where one is recorded), so that experiments recorded
with different settings can be compared. The noise of each plane is NMRTools'
whole-spectrum estimate (`spec[:noise]`), or where that is missing a robust estimate from
the plane itself.
"""
struct SpectrumSeries
    shifts::Vector{Float64}
    intensities::Matrix{Float64}
    noise::Vector{Float64}
    sources::Vector{String}
    experiments::Vector{Int}
    nucleus::String
    title::String
end

function SpectrumSeries(shifts, intensities, noise, sources, experiments; nucleus="1H",
                        title="")
    n = size(intensities, 2)
    length(shifts) == size(intensities, 1) ||
        throw(DimensionMismatch("$(length(shifts)) shifts for $(size(intensities, 1)) points"))
    length(noise) == length(sources) == length(experiments) == n ||
        throw(DimensionMismatch("noise, sources and experiments must have one entry per plane"))
    return SpectrumSeries(collect(Float64, shifts), Matrix{Float64}(intensities),
                          collect(Float64, noise), collect(String, sources),
                          collect(Int, experiments), String(nucleus), String(title))
end

function SpectrumSeries(spectra::AbstractVector{<:NMRData})
    isempty(spectra) && throw(ArgumentError("no spectra to track peaks in"))
    blocks = planesof.(spectra)
    δ = first(blocks).shifts
    for (k, b) in enumerate(blocks)
        sameaxis(δ, b.shifts) ||
            throw(DimensionMismatch("the chemical shift axis of spectrum $k differs from the " *
                                    "first; every spectrum must share one axis"))
    end
    experiments = reduce(vcat, [fill(k, size(b.y, 2)) for (k, b) in enumerate(blocks)])
    spec = frequencyfirst(first(spectra))
    return SpectrumSeries(δ, reduce(hcat, [b.y for b in blocks]),
                          reduce(vcat, [b.noise for b in blocks]),
                          reduce(vcat, [fill(b.source, size(b.y, 2)) for b in blocks]),
                          experiments;
                          nucleus=string(something(label(dims(spec, 1)), "1H")),
                          title=string(something(label(spec), "")))
end

"""
The planes of one processed spectrum: its shift axis, real intensities and noise, both
divided by `scale(spec)`, and its source.
"""
function planesof(spec::NMRData)
    spec = frequencyfirst(spec)
    δ = collect(Float64, data(dims(spec, 1)))
    f = scale(spec)
    y = Matrix{Float64}(reshape(real.(data(spec)), length(δ), :)) ./ f
    σ = spec[:noise]
    noise = isnothing(σ) ? map(noiselevel, eachcol(y)) : fill(Float64(σ / f), size(y, 2))
    return (; shifts=δ, y, noise, source=string(something(spec[:filename], spec[:title], "")))
end

frequencyfirst(spec) = frequencyfirst(spec, dims(spec, 1))
frequencyfirst(spec, ::FrequencyDimension) = spec
frequencyfirst(spec, ::Any) = permutedims(spec, (2, 1))

"Do two shift axes agree, to within a thousandth of a point?"
function sameaxis(a, b)
    length(a) == length(b) || return false
    tol = 1e-3 * abs(a[2] - a[1])
    return abs(first(a) - first(b)) ≤ tol && abs(last(a) - last(b)) ≤ tol
end

"Robust rms noise estimate: median absolute deviation scaled for a Gaussian."
noiselevel(y) = 1.4826 * median(abs.(y .- median(y)))

Base.length(s::SpectrumSeries) = size(s.intensities, 2)
Base.eachindex(s::SpectrumSeries) = axes(s.intensities, 2)
shifts(s::SpectrumSeries) = s.shifts
intensities(s::SpectrumSeries, i) = view(s.intensities, :, i)
noiselevel(s::SpectrumSeries, i) = s.noise[i]
sources(s::SpectrumSeries) = s.sources
"The planes processed from experiment `k`."
planesfrom(s::SpectrumSeries, k) = findall(==(k), s.experiments)

"""
    displayscale(s; snr=5) -> Float64

The intensity that fills one row of the stack at unit gain: the 90th percentile of the
heights of the peaks above `snr` times the noise, across every plane. A percentile of peak
heights rather than the largest point, so that residual water or a concentrated buffer does
not flatten everything else when the whole spectrum is shown.
"""
function displayscale(s::SpectrumSeries; snr=5)
    heights = Float64[]
    for i in eachindex(s)
        y = intensities(s, i)
        threshold = snr * noiselevel(s, i)
        for k in 2:(length(y) - 1)
            y[k] ≥ y[k - 1] && y[k] > y[k + 1] && y[k] ≥ threshold && push!(heights, y[k])
        end
    end
    isempty(heights) && return max(maximum(abs, s.intensities), eps())
    return quantile(heights, 0.9)
end
