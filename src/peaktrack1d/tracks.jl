# The picking model: multiplets made of lines, and lines made of picks, at most one per
# plane. No Makie and no NMRData.

"""
    Pick(plane, shift, height)

A line's position in one plane: the chemical shift as measured (ppm, never aligned) and the
intensity there.
"""
struct Pick
    plane::Int
    shift::Float64
    height::Float64
end

"""
    Multiplet(label, lines)

A named multiplet and the lines it is made of, by line number. A singlet is a multiplet of
one line. Lines are referred to by number so that renaming or regrouping never touches
their picks.
"""
struct Multiplet
    label::String
    lines::Vector{Int}
end

"""
    Tracks(labels)

Every multiplet, the reference first, and every pick, keyed by `(line, plane)`. Starts with
one multiplet of one empty line per label.
"""
mutable struct Tracks
    multiplets::Vector{Multiplet}
    picks::Dict{Tuple{Int,Int},Pick}
    nextline::Int
end

function Tracks(labels::AbstractVector{<:AbstractString})
    isempty(labels) && throw(ArgumentError("at least a reference label is needed"))
    multiplets = [Multiplet(String(l), [k]) for (k, l) in enumerate(labels)]
    return Tracks(multiplets, Dict{Tuple{Int,Int},Pick}(), length(labels) + 1)
end

function Base.copy(t::Tracks)
    return Tracks([Multiplet(m.label, copy(m.lines)) for m in t.multiplets], copy(t.picks),
                  t.nextline)
end

"Index of the reference multiplet, shown to the user as 0."
const REF = 1

"The reference multiplet's name unless another is given."
const DEFAULTREFERENCE = "reference"

multiplets(t::Tracks) = t.multiplets
nmultiplets(t::Tracks) = length(t.multiplets)
multipletlabels(t::Tracks) = [m.label for m in t.multiplets]
linesof(t::Tracks, m) = t.multiplets[m].lines

"Index of the multiplet holding line `l`, or `nothing`."
multipletof(t::Tracks, l) = findfirst(m -> l in m.lines, t.multiplets)
hasline(t::Tracks, l) = !isnothing(multipletof(t, l))
isreference(t::Tracks, l) = multipletof(t, l) == REF

pickat(t::Tracks, l, i) = get(t.picks, (l, i), nothing)
picksof(t::Tracks, l) = sort!([q for ((k, _), q) in t.picks if k == l]; by=q -> q.plane)
haspicks(t::Tracks, l) = any(k -> first(k) == l, keys(t.picks))
haspicks(t::Tracks) = !isempty(t.picks)
npicks(t::Tracks) = length(t.picks)

"Mean measured shift of a line over the planes where it is picked; `NaN` for an empty line."
function meanshift(t::Tracks, l)
    qs = picksof(t, l)
    return isempty(qs) ? NaN : mean(q.shift for q in qs)
end

"""
    orderedlines(t, m) -> Vector{Int}

The lines of multiplet `m` that have been picked anywhere, in decreasing order of shift, as
a spectrum reads from left to right. Their position in this list is the line number written
to the output files.
"""
function orderedlines(t::Tracks, m)
    picked = filter(l -> haspicks(t, l), linesof(t, m))
    return sort!(picked; by=l -> meanshift(t, l), rev=true)
end

"Every line of multiplet `m` in the order they are stepped through: picked ones first."
navorder(t::Tracks, m) = [orderedlines(t, m); filter(l -> !haspicks(t, l), linesof(t, m))]

"""
    centre(t, m, i) -> Union{Float64,Nothing}

The chemical shift of multiplet `m` in plane `i`: the mean of its lines' shifts, which is
the centre of any symmetric first-order multiplet. Defined only where every line picked
anywhere is picked in this plane, since a mean over whichever lines happen to be present
would be biased in exactly the planes where one is missing.
"""
centre(t::Tracks, m, i) = linecentre(t, orderedlines(t, m), i)

function linecentre(t::Tracks, lines, i)
    isempty(lines) && return nothing
    qs = [pickat(t, l, i) for l in lines]
    any(isnothing, qs) && return nothing
    return mean(q.shift for q in qs)
end

"The reference shift in each of `n` planes, or `nothing` if any plane lacks it."
function referenceshifts(t::Tracks, n)
    δ = [centre(t, REF, i) for i in 1:n]
    return any(isnothing, δ) ? nothing : Float64.(δ)
end

"""
    vacantname(labels) -> String

The first of `peak1`, `peak2`, ... not already in `labels`, the name a multiplet takes when
it is created by picking.
"""
function vacantname(labels)
    return first(n for n in ("peak$k" for k in Iterators.countfrom(1)) if n ∉ labels)
end

# ---------------------------------------------------------------- editing

"Add a multiplet of one new, empty line called `name`, returning the line's number."
function addmultiplet!(t::Tracks, name)
    l = t.nextline
    t.nextline += 1
    push!(t.multiplets, Multiplet(String(name), [l]))
    return l
end

"Place `q` on line `l`, replacing any earlier pick of that line in the same plane."
setpick!(t::Tracks, l, q::Pick) = (t.picks[(l, q.plane)]=q; t)
setpick!(t::Tracks, l, ::Nothing) = t

"Remove the picks of line `l` for which `f(pick)` is true."
function erase!(f, t::Tracks, l)
    filter!(kv -> first(kv.first) != l || !f(kv.second), t.picks)
    return t
end

"""
    group!(t, a, b) -> Int

Join the multiplet holding line `b` to the one holding line `a`, keeping the label of `a`'s
multiplet, unless `b` belongs to the reference, which is never absorbed into another
multiplet. Returns the index of the joined multiplet.
"""
function group!(t::Tracks, a, b)
    ma, mb = multipletof(t, a), multipletof(t, b)
    (isnothing(ma) || isnothing(mb)) && throw(ArgumentError("no line $a or $b to group"))
    ma == mb && return ma
    keep, drop = mb == REF ? (mb, ma) : (ma, mb)
    append!(t.multiplets[keep].lines, t.multiplets[drop].lines)
    deleteat!(t.multiplets, drop)
    return drop < keep ? keep - 1 : keep
end

"""
    ungroup!(t, l) -> Int

Move line `l` out of its multiplet into a new multiplet of its own, placed just after it.
Returns the index of the multiplet that now holds `l`.
"""
function ungroup!(t::Tracks, l)
    m = multipletof(t, l)
    isnothing(m) && throw(ArgumentError("no line $l to ungroup"))
    length(linesof(t, m)) == 1 && return m
    filter!(!=(l), t.multiplets[m].lines)
    insert!(t.multiplets, m + 1, Multiplet(vacantname(multipletlabels(t)), [l]))
    return m + 1
end

"""
    deleteline!(t, l)

Remove line `l` and its picks. A multiplet left without lines goes too, except the
reference, which keeps one empty line so that it can be picked again.
"""
function deleteline!(t::Tracks, l)
    m = multipletof(t, l)
    isnothing(m) && return t
    erase!(_ -> true, t, l)
    if length(linesof(t, m)) > 1
        filter!(!=(l), t.multiplets[m].lines)
    elseif m != REF
        deleteat!(t.multiplets, m)
    end
    return t
end

"""
    cleanlabel(name) -> String

A multiplet label as it is stored: trimmed, without commas, since labels are written
unquoted to CSV, and with underscores replaced by hyphens, since in a CSV header an
underscore separates a quantity from a tag (see `docs/src/advanced/conventions.md`).
"""
cleanlabel(name::AbstractString) = replace(sanitizelabel(strip(name)), '_' => '-')

"""
    relabel!(t, m, name) -> Bool

Rename multiplet `m`, cleaning the name with [`cleanlabel`](@ref). Returns `false`,
changing nothing, if the name is blank or already taken by another multiplet.
"""
function relabel!(t::Tracks, m, name)
    name = cleanlabel(name)
    isempty(name) && return false
    name != t.multiplets[m].label && name in multipletlabels(t) && return false
    t.multiplets[m] = Multiplet(name, t.multiplets[m].lines)
    return true
end
