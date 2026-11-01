# The analysis a window returns, and the tables made from it.

"""
    PeakTracker

A peak-tracking analysis, returned by [`peaktrack1d`](@ref) when its window closes: the
spectra, the coordinates distinguishing the planes, the multiplets and their picks, how the
spectra were processed (for time-domain data), and how the shifts are referenced (see
[`Referencing`](@ref)), along with the window it was picked in.

Show it in the REPL for a table of multiplet centres, copy every line's shifts with
`clipboard(p)`, get every pick with `planeresults(p)`, draw the stack with `stackplot(p)`, or reopen the window with
`peaktrack1d(p)`.
"""
struct PeakTracker
    series::SpectrumSeries
    coordinates::NamedTuple
    tracks::Observable{Tracks}
    processing::Union{Nothing,Processing}
    reference::Referencing
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
function coordinatetable(p::Pair, n)
    return coordinatetable(NamedTuple{(Symbol(first(p)),)}((last(p),)), n)
end
function coordinatetable(nt::NamedTuple, n)
    for (name, v) in pairs(nt)
        isnothing(n) || length(v) == n ||
            throw(DimensionMismatch("$(length(v)) values of $name for $n spectra"))
    end
    return map(collect, nt)
end

"""
    referencedshift(p, i, δ) -> Union{Float64,Nothing}

The measured shift `δ` (ppm) in plane `i`, referenced as `p.reference` says, or `nothing`
where the reference is missing in that plane.
"""
referencedshift(p::PeakTracker, i, δ) = referencedshift(p.reference, p, i, δ)
referencedshift(p::PeakTracker, q::Pick) = referencedshift(p, q.plane, q.shift)

function referencedshift(r::PickedReference, p::PeakTracker, i, δ)
    c = centre(tracks(p), REF, i)
    return isnothing(c) ? nothing : δ - c + r.shift
end

function referencedshift(r::IndirectReference, p::PeakTracker, i, δ)
    ν0 = r.frequencies[i]
    isnan(ν0) && return nothing
    s = p.series
    return 1e6 * (s.bf[i] * (1 + 1e-6δ) / (s.xi * ν0) - 1)
end

"The centre of the multiplet made of `lines` in plane `i`, referenced; `nothing` where undefined."
function referenced(p::PeakTracker, lines, i)
    c = linecentre(tracks(p), lines, i)
    return isnothing(c) ? nothing : referencedshift(p, i, c)
end

"Is the reference available in any plane, and in every plane?"
function isreferenced(p::PeakTracker)
    return any(i -> !isnothing(referencedshift(p, i, 0.0)), 1:nplanes(p))
end
function allreferenced(p::PeakTracker)
    return all(i -> !isnothing(referencedshift(p, i, 0.0)), 1:nplanes(p))
end

"""
    alignment(p) -> Union{Vector{Float64},Nothing}

How far to displace each plane of the stack so that referenced shifts line up, or `nothing`
unless the reference is available in every plane. Taken at the middle of the axis, which for
indirect referencing differs from anywhere else on it by parts in 10⁶ of the displacement.
"""
function alignment(p::PeakTracker)
    δ = sum(extrema(shifts(p.series))) / 2
    offsets = [referencedshift(p, i, δ) for i in 1:nplanes(p)]
    return any(isnothing, offsets) ? nothing : Float64.(offsets) .- δ
end

"""
    IndirectReference(source::PeakTracker, s::SpectrumSeries)

Reference the spectra `s` indirectly to the reference picked in the ¹H analysis `source`,
plane by plane. The two must have the same number of planes, matched in order, recorded at
the same field, and `source` must have its own picked reference, normally DSS at 0 ppm.
"""
function IndirectReference(source::PeakTracker, s::SpectrumSeries)
    h = source.series
    source.reference isa PickedReference ||
        throw(ArgumentError("the analysis to reference from must have a reference picked " *
                            "in it, not one referenced indirectly itself"))
    h.xi == 1 ||
        throw(ArgumentError("indirect referencing needs a ¹H analysis with DSS picked, " *
                            "not a $(h.nucleus) one"))
    isfinite(s.xi) ||
        throw(ArgumentError("NMRTools has no Ξ ratio for $(s.nucleus), so it can't be " *
                            "referenced indirectly"))
    length(h) == length(s) ||
        throw(DimensionMismatch("the ¹H analysis has $(length(h)) spectra and this one " *
                                "$(length(s)); they must match one to one, in order"))
    for i in eachindex(s)
        isfinite(h.bf[i]) && isfinite(s.bf[i]) ||
            throw(ArgumentError("the base frequency of spectrum $i is not known"))
        # the field is judged from the acquisition base frequencies, which the processing
        # reference SR does not move, wherever both are recorded
        fx, fh = isfinite(s.bf1[i]) && isfinite(h.bf1[i]) ? (s.bf1[i], h.bf1[i]) :
                 (s.bf[i], h.bf[i])
        abs(fx / (s.xi * fh) - 1) < FIELDTOLERANCE ||
            throw(ArgumentError("spectrum $i was recorded at a different field from ¹H " *
                                "spectrum $i: $(round(1e-6fx; digits=3)) MHz " *
                                "$(s.nucleus) against $(round(1e-6fh; digits=3)) MHz ¹H"))
    end
    t = tracks(source)
    r0 = source.reference.shift
    frequencies = map(eachindex(h)) do i
        c = centre(t, REF, i)
        return isnothing(c) ? NaN : h.bf[i] * (1 + 1e-6c) / (1 + 1e-6r0)
    end
    missingspectra = findall(isnan, frequencies)
    isempty(missingspectra) ||
        @warn "The ¹H reference is not picked in every spectrum, so these are not " *
              "referenced" missingspectra
    return IndirectReference(frequencies, multiplets(t)[REF].label)
end

"""
    referencing(reference, refshift, s) -> (label, Referencing)

The reference multiplet's label and how shifts are referenced, from the `reference`
keyword: a name picks the reference in this analysis, placed at `refshift`; a ¹H
`PeakTracker` references indirectly to the reference picked in it.
"""
referencing(name::AbstractString, refshift, s) = cleanlabel(name), PickedReference(refshift)
function referencing(source::PeakTracker, refshift, s)
    r = IndirectReference(source, s)
    return r.label, r
end
referencing(r::Tuple{AbstractString,Referencing}, refshift, s) = String(first(r)), last(r)

"How the shifts are referenced, in words."
function referencetext(p::PeakTracker)
    label = multiplets(tracks(p))[REF].label
    p.reference isa PickedReference &&
        return "Referenced to $label, placed at $(p.reference.shift) ppm"
    return "Referenced indirectly to $label in a 1H analysis of the same samples, on the " *
           "IUPAC unified scale (Xi = $(p.series.xi))"
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
    append!(header,
            [csvcolumn("shift", "ppm"), "height", "snr",
             csvcolumn("centre", "ppm"), csvcolumn("referenced", "ppm")])
    rows = Vector{String}[]
    for (m, mp) in enumerate(multiplets(t))
        ls = orderedlines(t, m)
        for (j, l) in enumerate(ls), i in eachindex(s)
            q = pickat(t, l, i)
            row = [cellvalue(sources(s)[i]), mp.label, string(j), string(i)]
            append!(row, [cellvalue(v[i]) for v in values(p.coordinates)])
            append!(row,
                    isnothing(q) ? fill("NA", 3) :
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

# ---------------------------------------------------------------- summary panel

"The first numeric coordinate and its name, or the spectrum numbers and `\"spectrum\"`."
function xcoordinate(p::PeakTracker)
    name = findfirst(v -> eltype(v) <: Real, p.coordinates)
    isnothing(name) && return collect(1.0:nplanes(p)), "spectrum"
    return Float64.(p.coordinates[name]), string(name)
end

"""
    summaryseries(p, line) -> NamedTuple

What the summary panel shows for the multiplet holding `line`, against the first numeric
coordinate or the plane number: the shift of each line and the multiplet's centre, and the
height of each line. Shifts are referenced where the reference is available in every plane,
and as measured otherwise, so one panel never mixes the two. Points run in order of the
coordinate, with `NaN` where a plane has no value and between lines, so the drawn lines
break there. The vacant slot (`line` 0) gives empty series and `m == 0`.
"""
function summaryseries(p::PeakTracker, line)
    t, n = tracks(p), nplanes(p)
    x, xlabel = xcoordinate(p)
    isref = allreferenced(p)
    name = referencename(t)
    ylabel = !isref ? "shift as measured (ppm)" :
             name == "the reference peak" ? "referenced shift (ppm)" :
             "shift referenced to $name (ppm)"
    m = something(multipletof(t, line), 0)
    if m == 0
        return (; xlabel, ylabel, m, label="", lines=Point2f[], centre=Point2f[],
                heights=Point2f[])
    end
    order = sortperm(x)
    ls = orderedlines(t, m)
    shift(δ, i) = isref ? referencedshift(p, i, δ) : δ
    gap = Point2f(NaN, NaN)
    function perline(value)
        # push! rather than [points; gap]: a Point2f is itself a vector, which vcat would
        # splice in as two numbers
        return reduce(vcat,
                      [push!([Point2f(x[i], value(pickat(t, l, i), i)) for i in order],
                             gap)
                       for l in ls]; init=Point2f[])
    end
    lines = perline((q, i) -> isnothing(q) ? NaN : shift(q.shift, i))
    heights = perline((q, i) -> isnothing(q) ? NaN : q.height)
    centre = map(order) do i
        c = linecentre(t, ls, i)
        return Point2f(x[i], isnothing(c) ? NaN : shift(c, i))
    end
    return (; xlabel, ylabel, m, label=multiplets(t)[m].label, lines, centre, heights)
end

# ---------------------------------------------------------------- display and clipboard

"""
The reference as a sentence names it: by its own name, such as DSS, or as "the reference
peak" while it keeps the default name.
"""
function referencename(t::Tracks)
    label = multiplets(t)[REF].label
    return label == DEFAULTREFERENCE ? "the reference peak" : label
end

"Spectra `ks` of `n`, as a sentence lists them: by number when there are a few, as a count otherwise."
function spectralist(ks, n)
    length(ks) == 1 && return "spectrum $(only(ks))"
    length(ks) ≤ 6 && return "spectra $(join(ks, ", ", " and "))"
    return "$(length(ks)) of the $n spectra"
end

"A shift as a sentence gives it: `0 ppm` rather than `0.0 ppm`."
ppmtext(δ) = isinteger(δ) ? "$(Int(δ)) ppm" : "$δ ppm"

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
    what = isreferenced(p) ? "referenced to $(multiplets(tracks(p))[REF].label)" :
           "as measured"
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

roundcell(::Nothing; digits) = ""
function roundcell(x::AbstractFloat; digits)
    return isinteger(x) ? string(Int(x)) : string(round(x; digits))
end
roundcell(x; digits) = string(x)

"""
    widetable(p; digits=4) -> (header, cells)

Multiplet centres as measured, one row per plane and one column per multiplet, the
reference first: the table shown in the REPL, rounded and with blank cells where a centre
is undefined, for reading rather than computing with.
"""
function widetable(p::PeakTracker; digits=4)
    t, n = tracks(p), nplanes(p)
    header = String["spectrum"; [string(k) for k in keys(p.coordinates)]]
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
