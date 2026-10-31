# The analysis a window returns, and the tables made from it.

"""
    PeakTracker

A peak-tracking analysis, returned by [`peaktrack1d`](@ref) when its window closes: the
spectra, the coordinates distinguishing the planes, the multiplets and their picks, how the
spectra were processed (for time-domain data), and the reference settings, along with the
window it was picked in.

Show it in the REPL for a table of multiplet centres, copy that table with `clipboard(p)`,
get every pick with `planeresults(p)`, draw the stack with `stackplot(p)`, or reopen the
window with `peaktrack1d(p)`.
"""
struct PeakTracker
    series::SpectrumSeries
    coordinates::NamedTuple
    tracks::Observable{Tracks}
    processing::Union{Nothing,Processing}
    refshift::Float64
    aligned::Observable{Bool}
    view::StackView
    figure::Figure
    axis::Axis
    call::Union{Nothing,AnalysisCall}
    outputdir::Observable{String}
end

tracks(p::PeakTracker) = p.tracks[]
nplanes(p::PeakTracker) = length(p.series)

"""
    coordinatetable(coordinates, n) -> NamedTuple

The coordinates of `n` planes as a `NamedTuple` of vectors, from a `NamedTuple` of vectors,
a `name => values` pair, a plain vector (named `x`), or `nothing` (none). With `n` nothing,
the lengths are not checked.
"""
coordinatetable(::Nothing, n) = NamedTuple()
coordinatetable(v::AbstractVector, n) = coordinatetable((; x=v), n)
coordinatetable(p::Pair, n) = coordinatetable(NamedTuple{(Symbol(first(p)),)}((last(p),)), n)
function coordinatetable(nt::NamedTuple, n)
    for (name, v) in pairs(nt)
        isnothing(n) || length(v) == n ||
            throw(DimensionMismatch("$(length(v)) values of $name for $n planes"))
    end
    return map(collect, nt)
end

"The reference shift in plane `i` (ppm), or `nothing`."
refcentre(p::PeakTracker, i) = centre(tracks(p), REF, i)

"""
The centre of the multiplet made of `lines` in plane `i`, referenced to the reference
multiplet placed at `refshift`; `nothing` where either is undefined.
"""
function referenced(p::PeakTracker, lines, i)
    c, r = linecentre(tracks(p), lines, i), refcentre(p, i)
    (isnothing(c) || isnothing(r)) && return nothing
    return c - r + p.refshift
end

# ---------------------------------------------------------------- long tables

cellvalue(x::AbstractString) = sanitizelabel(x)
cellvalue(x) = csvvalue(x)

"""
    seriestable(p) -> (header, rows)

Column names and rows of `series.csv`: one row per line per plane, keyed by multiplet
`label`, `line` (numbered by decreasing shift) and `plane`, with the plane's coordinates,
the line's measured `shift` (ppm), `height` and `snr`, and the multiplet's `centre` and
`referenced` shift (ppm) in that plane, repeated down its lines. Absent values are `NA`.
"""
function seriestable(p::PeakTracker)
    t, s = tracks(p), p.series
    header = ["source", "label", "line", "plane"]
    append!(header, string.(keys(p.coordinates)))
    append!(header, [csvcolumn("shift", "ppm"), "height", "snr",
                     csvcolumn("centre", "ppm"), csvcolumn("referenced", "ppm")])
    rows = Vector{String}[]
    for (m, mp) in enumerate(multiplets(t))
        ls = orderedlines(t, m)
        for (j, l) in enumerate(ls), i in eachindex(s)
            q = pickat(t, l, i)
            row = [cellvalue(sources(s)[i]), mp.label, string(j), string(i)]
            append!(row, [cellvalue(v[i]) for v in values(p.coordinates)])
            append!(row, isnothing(q) ? fill("NA", 3) :
                         [csvvalue(q.shift), csvvalue(q.height), csvvalue(snr(s, q))])
            push!(row, csvvalue(linecentre(t, ls, i)), csvvalue(referenced(p, ls, i)))
            push!(rows, row)
        end
    end
    return header, rows
end

"""
    peaklisttable(p) -> (header, rows)

Column names and rows of `peaklist.csv`, what the user placed: one row per pick, keyed by
multiplet `label`, whether it is the `reference`, `line` and `plane`. A multiplet with no
picks has one row with blank `line` and `plane` and an `NA` shift, so that its name survives
a reload.
"""
function peaklisttable(p::PeakTracker)
    t = tracks(p)
    header = ["label", "reference", "line", "plane", csvcolumn("shift", "ppm")]
    rows = Vector{String}[]
    for (m, mp) in enumerate(multiplets(t))
        isref = string(m == REF)
        ls = orderedlines(t, m)
        isempty(ls) && push!(rows, [mp.label, isref, "", "", "NA"])
        for (j, l) in enumerate(ls), q in picksof(t, l)
            push!(rows, [mp.label, isref, string(j), string(q.plane), csvvalue(q.shift)])
        end
    end
    return header, rows
end

"""
    planeresults(p::PeakTracker) -> Vector{NamedTuple}

One row per line per plane, the counterpart of `series.csv`: the multiplet `label`, `line`
and `plane`, the plane's coordinates, the line's `shift` (ppm), `height` and `snr`, and the
multiplet's `centre` and `referenced` shift (ppm), with `missing` where absent.
"""
function planeresults(p::PeakTracker)
    t, s = tracks(p), p.series
    orm(x) = something(x, missing)
    rows = NamedTuple[]
    for (m, mp) in enumerate(multiplets(t))
        ls = orderedlines(t, m)
        for (j, l) in enumerate(ls), i in eachindex(s)
            q = pickat(t, l, i)
            coords = map(v -> v[i], p.coordinates)
            push!(rows,
                  merge((; label=mp.label, line=j, plane=i), coords,
                        (; shift=isnothing(q) ? missing : q.shift,
                         height=isnothing(q) ? missing : q.height,
                         snr=isnothing(q) ? missing : snr(s, q),
                         centre=orm(linecentre(t, ls, i)),
                         referenced=orm(referenced(p, ls, i)))))
        end
    end
    return rows
end

# ---------------------------------------------------------------- the wide table

roundcell(::Nothing; digits) = ""
roundcell(x::AbstractFloat; digits) = isinteger(x) ? string(Int(x)) : string(round(x; digits))
roundcell(x; digits) = string(x)

"""
    widetable(p; digits=4) -> (header, cells)

Multiplet centres as measured, one row per plane and one column per multiplet, the
reference first: the table shown in the REPL and copied to the clipboard, rounded and with
blank cells where a centre is undefined, for reading and pasting into a spreadsheet rather
than for computing with.
"""
function widetable(p::PeakTracker; digits=4)
    t, n = tracks(p), nplanes(p)
    header = ["plane"; string.(keys(p.coordinates))]
    columns = Vector{String}[string.(1:n)]
    for v in values(p.coordinates)
        push!(columns, [roundcell(x; digits) for x in v])
    end
    for (m, mp) in enumerate(multiplets(t))
        push!(header, m == REF ? "$(mp.label) (reference)" : mp.label)
        push!(columns, [roundcell(centre(t, m, i); digits) for i in 1:n])
    end
    return header, reduce(hcat, columns)
end

function delimited(io::IO, p::PeakTracker, delim)
    header, cells = widetable(p)
    println(io, join(header, delim))
    for r in eachrow(cells)
        println(io, join(r, delim))
    end
end

Base.show(io::IO, ::MIME"text/tab-separated-values", p::PeakTracker) = delimited(io, p, '\t')

function Base.show(io::IO, p::PeakTracker)
    t = tracks(p)
    names = map(eachindex(multiplets(t))) do m
        k = length(orderedlines(t, m))
        return k > 1 ? "$(multiplets(t)[m].label) ($k lines)" : multiplets(t)[m].label
    end
    return print(io, "PeakTracker(", join(names, ", "), "; $(npicks(t)) picks)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PeakTracker)
    show(io, p)
    haspicks(tracks(p)) || return
    println(io)
    header, cells = widetable(p)
    w = [maximum(textwidth, [header[j]; cells[:, j]]) for j in eachindex(header)]
    print(io, join(lpad.(header, w), "  "))
    for r in eachrow(cells)
        print(io, "\n", join(lpad.(r, w), "  "))
    end
end

"Copy the table of multiplet centres to the clipboard, tab-separated, to paste into a spreadsheet."
function InteractiveUtils.clipboard(p::PeakTracker)
    return clipboard(sprint(show, MIME"text/tab-separated-values"(), p))
end
