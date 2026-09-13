"""
    Trace(δ, y)

A single 1D spectrum: a chemical-shift axis `δ` (ppm) and intensities `y`.

Plain vectors, with no dependency on NMRData, Makie or GUI state: the analysis layer
operates on `Trace`s, and the adapters in `nmrdata.jl` convert NMRData into them.
"""
struct Trace
    δ::Vector{Float64}
    y::Vector{Float64}
    function Trace(δ, y)
        length(δ) == length(y) ||
            throw(ArgumentError("δ and y must have equal length ($(length(δ)) vs $(length(y)))"))
        return new(collect(float.(δ)), collect(float.(y)))
    end
end

Base.length(t::Trace) = length(t.y)

"""
    Region(label, lo, hi)
    Region(label, δ)

A named chemical-shift interval (ppm) used for integration. `lo`/`hi` are stored sorted,
so the order in which the bounds are supplied does not matter. A zero-width region
(`lo == hi`, or the single-shift constructor) selects the nearest point — i.e. a peak
height.
"""
struct Region
    label::String
    lo::Float64
    hi::Float64
    function Region(label, a, b)
        return new(String(label), min(float(a), float(b)), max(float(a), float(b)))
    end
end

Region(label, δ::Real) = Region(label, δ, δ)

"""
    defaultregionwidth(δ) -> Float64

A default region width: 2% of the chemical-shift range spanned by `δ`, matching the noise
marker's drag-handle width.
"""
defaultregionwidth(δ::AbstractVector) = 0.02 * (maximum(δ) - minimum(δ))

width(r::Region) = r.hi - r.lo

"""
    centre(region) -> Float64

Midpoint of a region (ppm). Regions are stored by their bounds, but centre-and-width is
the natural handle for interaction, so both views are available.
"""
centre(r::Region) = (r.lo + r.hi) / 2

"""
    recentre(region, c) -> Region

The same region, same width and label, moved to be centred on `c`.
"""
recentre(r::Region, c) = Region(r.label, c - width(r) / 2, c + width(r) / 2)

"""
    setwidth(region, w) -> Region

The same region, same label and centre, with width `w`.
"""
setwidth(r::Region, w) = Region(r.label, centre(r) - w / 2, centre(r) + w / 2)

"""
    Planes(traces, vars)

Long-format collection of 1D spectra: one `Trace` per row, with `vars[i]` a `NamedTuple`
giving the values of the arrayed acquisition variables for that spectrum
(e.g. `(; time = 0.1, which = :trosy)`). All rows must share the same set of variable
names.
"""
struct Planes
    traces::Vector{Trace}
    vars::Vector{<:NamedTuple}
    function Planes(traces, vars)
        length(traces) == length(vars) ||
            throw(ArgumentError("traces and vars must have equal length ($(length(traces)) vs $(length(vars)))"))
        isempty(vars) || allequal(keys.(vars)) ||
            throw(ArgumentError("all planes must share the same variable names"))
        return new(collect(Trace, traces), collect(vars))
    end
end

nplanes(p::Planes) = length(p.traces)

"""
    column(planes, name) -> Vector

Return the values of variable `name` across all planes.
"""
column(p::Planes, name::Symbol) = [v[name] for v in p.vars]

"""
    hasvar(planes, name) -> Bool

Whether the planes carry an arrayed variable called `name`.
"""
hasvar(p::Planes, name::Symbol) = !isempty(p.vars) && haskey(first(p.vars), name)

"""
    Dataset1D(planes, noisecenter, [label, [sources]])

The planes plus the universal noise position (ppm), a `label` saying where the data came
from, and `sources`, the same per plane. Both are plain `String`s, so results files can
name what they were computed from without NMRData re-entering the analysis core.

`sources` is per plane because an experiment may combine files: TRACT loads a TROSY and an
anti-TROSY spectrum. It defaults to `label` repeated. It is not a plane variable, being
provenance rather than a coordinate: two files may be replicates with identical
coordinates.

Only the noise *centre* is stored. The noise region takes the width of whichever signal
region is being measured, which is what makes its integral a direct estimate of that
region's noise (see [`integrate`](@ref)).
"""
struct Dataset1D
    planes::Planes
    noisecenter::Float64
    label::String
    sources::Vector{String}
    function Dataset1D(planes, noisecenter, label, sources)
        length(sources) == nplanes(planes) ||
            throw(ArgumentError("got $(length(sources)) sources for " *
                                "$(nplanes(planes)) planes"))
        return new(planes, Float64(noisecenter), String(label), collect(String, sources))
    end
end

function Dataset1D(planes, noisecenter, label="")
    return Dataset1D(planes, noisecenter, label, fill(String(label), nplanes(planes)))
end

nplanes(d::Dataset1D) = nplanes(d.planes)

"""
    sources(dataset) -> Vector{String}

Where each plane came from, one entry per plane. This is what the `source` column of
`series.csv` carries.
"""
sources(d::Dataset1D) = d.sources

"""
    groupseries(planes, cols) -> Vector{Pair{NamedTuple,Vector{Int}}}

Group plane indices by the values of the variables named in `cols`, preserving order of
first appearance. With `cols = ()` every plane falls in one group keyed by the empty
`NamedTuple`. A series is the set of planes sharing all grouping variables, differing only
in the fit axis.
"""
function groupseries(planes::Planes, cols::Tuple)
    keyfor(i) = NamedTuple{cols}(Tuple(planes.vars[i][c] for c in cols))
    keys = [keyfor(i) for i in 1:nplanes(planes)]
    uniquekeys = unique(keys)
    return [k => findall(==(k), keys) for k in uniquekeys]
end
