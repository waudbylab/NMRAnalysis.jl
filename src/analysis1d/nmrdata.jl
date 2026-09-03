# NMRData → analysis-core adapters. These are the only functions in the module that touch
# NMRData; everything downstream works on plain vectors (`Trace`/`Planes`/`Dataset1D`).
# The per-experiment entry points that use them live in the `expt-*.jl` files.

loadspec(x::AbstractString) = loadnmr(String(x))
loadspec(x::Integer) = loadnmr(string(x))
loadspec(x) = x

"""Annotation lookup returning `nothing` rather than throwing when absent."""
function annotation(spec, keys...)
    try
        return annotations(spec, keys...)
    catch
        return nothing
    end
end

"""
    tracesfromspec(spec) -> Vector{Trace}

Extract one `Trace` per plane from an N-dimensional NMRData whose first axis is the
chemical shift (`F1Dim`) — a pseudo-2D (shift × planes), or a pseudo-3D/-nD such as an
`Exchange1D` off-resonance R1ρ spectrum (shift × spinlock × delay, say). The non-shift
axes are flattened (in the array's natural, column-major order) into one list of planes;
`Planes`/`Dataset1D` only ever need "one 1D trace per plane" regardless of how many
arrayed dimensions produced it.

Intensities are divided by the spectrum's estimated RMS noise level (`spec[:noise]`,
computed once by NMRTools when the spectrum is loaded) here, at the point traces are
built - not live, per-analysis, from a user-positioned noise region - so every
downstream `Trace` is already expressed in signal-to-noise units.
"""
function tracesfromspec(spec)
    δ = collect(data(spec, F1Dim))
    Y = data(spec)
    n = length(δ)
    size(Y, 1) == n ||
        throw(ArgumentError("expected the chemical shift (F1Dim) to be the first " *
                            "dimension of $(typeof(spec)); got size $(size(Y))"))
    noise = something(spec[:noise], 1.0)
    Yflat = reshape(Y, n, :) ./ noise
    return [Trace(δ, collect(view(Yflat, :, i))) for i in 1:size(Yflat, 2)]
end

"""
    defaultnoisecentre(spec; frac=0.9) -> Float64

A fallback noise position: `frac` of the way across the chemical-shift axis. The GUI lets
the user reposition it.
"""
function defaultnoisecentre(spec; frac=0.9)
    δ = collect(data(spec, F1Dim))
    lo, hi = extrema(δ)
    return lo + frac * (hi - lo)
end

"""
    datasetfromspec(spec, vars; noisecenter=defaultnoisecentre(spec)) -> Dataset1D

Build a `Dataset1D` from a pseudo-2D `spec` and a vector of per-plane variable
`NamedTuple`s (`length(vars) == number of planes`).
"""
function datasetfromspec(spec, vars::AbstractVector{<:NamedTuple};
                           noisecenter::Real=defaultnoisecentre(spec))
    return Dataset1D(Planes(tracesfromspec(spec), collect(vars)), Float64(noisecenter),
                     speclabel(spec))
end

"""A short description of where a spectrum came from, for results-file provenance."""
speclabel(spec) = string(something(spec[:filename], spec[:title], ""))
