"""
    Trace(δ, y)

A single 1D spectrum: a chemical-shift axis `δ` (ppm) and intensities `y`.

Deliberately holds plain vectors and has no dependency on NMRData, Makie, or any GUI
state — the analysis layer operates on `Trace`s so the science stays independent of the
GUI. Loaders (`loaders.jl`) convert NMRData into `Trace`s.
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
    Region(label, a, b) = new(String(label), min(float(a), float(b)), max(float(a), float(b)))
end

Region(label, δ::Real) = Region(label, δ, δ)

width(r::Region) = r.hi - r.lo

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
    Dataset1D(planes, noise)

The planes plus the universal noise `Region`. Noise handling is shared across all
experiments rather than re-implemented per analysis.
"""
struct Dataset1D
    planes::Planes
    noise::Region
end

nplanes(d::Dataset1D) = nplanes(d.planes)

"""
    groupseries(planes, cols) -> Vector{Pair{NamedTuple,Vector{Int}}}

Group plane indices by the values of the variables named in `cols` (a tuple of
`Symbol`s), preserving order of first appearance. With `cols = ()` every plane falls in
a single group keyed by the empty `NamedTuple`. This is the generic form of R1ρ's
`onresseries`: a "series" is the set of planes sharing all grouping variables, i.e.
differing only in the fit-axis.
"""
function groupseries(planes::Planes, cols::Tuple)
    keyfor(i) = NamedTuple{cols}(Tuple(planes.vars[i][c] for c in cols))
    keys = [keyfor(i) for i in 1:nplanes(planes)]
    uniquekeys = unique(keys)
    return [k => findall(==(k), keys) for k in uniquekeys]
end
