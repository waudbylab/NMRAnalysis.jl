# The analysis a window returns, and the tables made from it.

"""
    PeakTracker

A peak-tracking analysis, returned by [`peaktrack1d`](@ref) when its window closes: the
spectra, the coordinates distinguishing the planes, the multiplets and their picks, how the
spectra were processed (for time-domain data), and the reference settings, along with the
window it was picked in.

Show it in the REPL for a summary, copy every line's shifts with `clipboard(p)`, get every
pick with `planeresults(p)`, draw the stack with `stackplot(p)`, or reopen the window with
`peaktrack1d(p)`.
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
`referenced` shift (ppm) in that plane, repeated down its lines. Absent values are `NA`. A
height is the intensity at the line's maximum, divided by `scale(spec)`: a height, not an
integral or an amplitude.
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

# ---------------------------------------------------------------- wide files

"""
    linetable(p, value, unit="") -> (header, rows)

A wide table for a spreadsheet: one row per plane, keyed by `plane` and the coordinates,
then one column per line of every multiplet that has been picked, the reference first. A
singlet's column is named by its label, and the lines of a multiplet `label[1]`,
`label[2]`, ... from the highest shift. Each cell is `value(q)` for the line's pick `q` in
that plane, at full precision, and `NA` where the line is not picked or `value` returns
`nothing`.
"""
function linetable(p::PeakTracker, value, unit="")
    t, n = tracks(p), nplanes(p)
    header = String["plane"; [string(k) for k in keys(p.coordinates)]]
    lines = Int[]
    for (m, mp) in enumerate(multiplets(t))
        ls = orderedlines(t, m)
        for (j, l) in enumerate(ls)
            push!(header, csvcolumn(length(ls) == 1 ? mp.label : "$(mp.label)[$j]", unit))
            push!(lines, l)
        end
    end
    rows = map(1:n) do i
        row = String[string(i)]
        append!(row, [cellvalue(v[i]) for v in values(p.coordinates)])
        for l in lines
            q = pickat(t, l, i)
            push!(row, csvvalue(isnothing(q) ? nothing : value(q)))
        end
        return row
    end
    return header, rows
end

"""
The shift of pick `q` referenced to the reference multiplet placed at `refshift`, or
`nothing` where the reference is not picked in that plane.
"""
function referencedshift(p::PeakTracker, q::Pick)
    r = refcentre(p, q.plane)
    return isnothing(r) ? nothing : q.shift - r + p.refshift
end

# ---------------------------------------------------------------- display and clipboard

"Has the reference been picked in any plane?"
isreferenced(p::PeakTracker) = any(i -> !isnothing(refcentre(p, i)), 1:nplanes(p))

"""
    clipboardtext(p) -> String

What `clipboard(p)` copies, tab-separated: the contents of `referencedshifts.csv` once the
reference has been picked in any plane, and of `shifts.csv` before then.
"""
function clipboardtext(p::PeakTracker)
    header, rows = isreferenced(p) ? linetable(p, q -> referencedshift(p, q), "ppm") :
                   linetable(p, q -> q.shift, "ppm")
    return join([join(header, '\t'); [join(r, '\t') for r in rows]], '\n') * '\n'
end

"""
    clipboard(p::PeakTracker)

Copy the shift of every line in every plane to the clipboard, tab-separated, to paste into a
spreadsheet: referenced to the reference multiplet once it has been picked, as in
`referencedshifts.csv`, and as measured before then, as in `shifts.csv`.
"""
function InteractiveUtils.clipboard(p::PeakTracker)
    clipboard(clipboardtext(p))
    what = isreferenced(p) ? "referenced to $(multiplets(tracks(p))[REF].label)" : "as measured"
    @info "Copied the shift of every line, $what, to the clipboard"
end

function Base.show(io::IO, p::PeakTracker)
    t = tracks(p)
    names = map(eachindex(multiplets(t))) do m
        k = length(orderedlines(t, m))
        return k > 1 ? "$(multiplets(t)[m].label) ($k lines)" : multiplets(t)[m].label
    end
    return print(io, "PeakTracker(", join(names, ", "), "; $(npicks(t)) picks)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PeakTracker)
    t, n = tracks(p), nplanes(p)
    print(io, "PeakTracker: $(npicks(t)) picks in $n planes")
    for (m, mp) in enumerate(multiplets(t))
        ls = orderedlines(t, m)
        k = count(i -> !isnothing(linecentre(t, ls, i)), 1:n)
        lines = length(ls) == 1 ? "1 line" : "$(length(ls)) lines"
        print(io, "\n  ", rpad(mp.label, 16), rpad(lines, 10), "centre in $k of $n planes",
              m == REF ? "  (reference)" : "")
    end
end
