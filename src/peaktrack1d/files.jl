# The output folder, and reading back what it holds. See docs/src/advanced/conventions.md.

"""
    readtable(path) -> (header, rows)

A CSV written by [`writetable`](@ref NMRAnalysis.writetable): `#` comment lines and blank
lines are skipped, the first remaining line is the header, and every field is stripped.
"""
function readtable(path::AbstractString)
    lines = filter(l -> !isempty(strip(l)) && !startswith(lstrip(l), '#'),
                   rstrip.(readlines(path), '\r'))
    isempty(lines) && throw(ArgumentError("$path holds no table"))
    fields(line) = String.(strip.(split(line, ',')))
    return fields(first(lines)), fields.(lines[2:end])
end

"Index of the column called `name`, or an error naming the file and the column."
function columnindex(header, name, path)
    j = findfirst(==(name), header)
    isnothing(j) && throw(ArgumentError("$path has no \"$name\" column"))
    return j
end

"""
    restore(s, header, rows; reference="reference", path="the peak list") -> Tracks

Multiplets and picks from the rows of a `peaklist.csv`. Each pick is placed exactly at its
saved shift, with its height read from `s`, so a peak list carries over to reprocessed or
different data. Lines are numbered within their multiplet as they were saved. The multiplet
marked as the reference comes first; without one, an empty reference called `reference`
is added.
"""
function restore(s::SpectrumSeries, header, rows; reference=DEFAULTREFERENCE,
                 path="the peak list")
    jl = columnindex(header, "label", path)
    jr = columnindex(header, "reference", path)
    jn = columnindex(header, "line", path)
    jp = columnindex(header, "plane", path)
    js = columnindex(header, csvcolumn("shift", "ppm"), path)

    labels = String[]
    reflabel = nothing
    for row in rows
        lbl = cleanlabel(row[jl])
        lbl in labels || push!(labels, lbl)
        lowercase(row[jr]) == "true" && isnothing(reflabel) && (reflabel = lbl)
    end
    if isnothing(reflabel)
        reflabel = cleanlabel(reference)
        reflabel in labels && throw(ArgumentError("$path marks no reference, and has a " *
                                                  "multiplet called $reflabel"))
    end
    t = Tracks([Multiplet(l, Int[]) for l in [reflabel; filter(!=(reflabel), labels)]],
               Dict{Tuple{Int,Int},Pick}(), 1)

    linenumbers = Dict{Tuple{String,String},Int}()
    for row in rows
        lbl = cleanlabel(row[jl])
        m = findfirst(==(lbl), multipletlabels(t))
        key = (lbl, row[jn])
        l = get!(linenumbers, key) do
            n = t.nextline
            t.nextline += 1
            push!(t.multiplets[m].lines, n)
            return n
        end
        i, δ = tryparse(Int, row[jp]), tryparse(Float64, row[js])
        (isnothing(i) || isnothing(δ)) && continue
        if !(1 ≤ i ≤ length(s))
            @warn "Spectrum $i of $lbl in $path is not in the data, so it was skipped"
            continue
        end
        setpick!(t, l, exact(s, i, δ))
    end
    for m in eachindex(t.multiplets)
        isempty(linesof(t, m)) || continue
        push!(t.multiplets[m].lines, t.nextline)
        t.nextline += 1
    end
    return t
end

function restore(s::SpectrumSeries, path::AbstractString; kw...)
    return restore(s, readtable(path)...; path, kw...)
end

"The `(p0, p1)` phases, in experiment order, from a `processing.csv`."
function readphases(path::AbstractString)
    header, rows = readtable(path)
    j0 = columnindex(header, csvcolumn("p0", "deg"), path)
    j1 = columnindex(header, csvcolumn("p1", "deg"), path)
    return [(parse(Float64, r[j0]), parse(Float64, r[j1])) for r in rows]
end

"""
    processingtable(p::Processing) -> (header, rows)

Column names and rows of `processing.csv`: one row per experiment, with where its phases
came from (`phasing`), the phases (degrees), and the line broadening, zero-filling and
solvent suppression. Passing the file back as `phase` repeats the processing without
phasing again.
"""
function processingtable(p::Processing)
    header = ["source", "phasing", csvcolumn("p0", "deg"), csvcolumn("p1", "deg"),
              csvcolumn("lb", "Hz"), "zerofill", "solvent"]
    rows = [[cellvalue(source(p, k)), string(p.phasing[k]), csvvalue(p.phases[k][1]),
             csvvalue(p.phases[k][2]), csvvalue(p.lb[]), string(p.zerofill),
             cellvalue(string(p.solvent))] for k in eachindex(p.fids)]
    return header, rows
end

"The comment lines heading every CSV: what was analysed, and from what."
function analysisinfo(p::PeakTracker)
    s = p.series
    lines = ["NMRAnalysis.jl $(pkgversion(NMRAnalysis))", "Analysis: peak tracking (1D)"]
    isempty(s.title) || push!(lines, "Title: $(s.title)")
    for src in unique(sources(s))
        isempty(src) || push!(lines, "Source: $src")
    end
    push!(lines, "Number of spectra: $(length(s))", "Nucleus: $(s.nucleus)")
    isempty(p.coordinates) ||
        push!(lines, "Coordinates: $(join(keys(p.coordinates), ", "))")
    push!(lines, referencetext(p))
    return lines
end

fmt(x; digits=4) = isnothing(x) ? "NA" : string(round(x; digits))

"""
    writesummary(path, p, folder)

Write `summary.txt`, the human-readable record, rounded for reading: what was analysed, how
it was processed, the reference, each multiplet's lines and the range of its centre, and
the call that repeats the analysis from the files saved in `folder`.
"""
function writesummary(path::AbstractString, p::PeakTracker, folder)
    t, n = tracks(p), nplanes(p)
    open(path, "w") do io
        foreach(line -> println(io, line), analysisinfo(p))
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        if !isnothing(p.processing)
            proc = p.processing
            println(io)
            println(io,
                    "Processing: exponential line broadening $(proc.lb[]) Hz, " *
                    "zero-filling ×$(proc.zerofill), solvent suppression $(proc.solvent)")
            println(io, "Phases (degrees, p0 at the first point):")
            for k in eachindex(proc.fids)
                p0, p1 = proc.phases[k]
                println(io,
                        "  $(rpad(shortpath(source(proc, k)), 30))  p0 $(lpad(fmt(p0; digits=1), 7))" *
                        "  p1 $(lpad(fmt(p1; digits=1), 7))  ($(proc.phasing[k]))")
            end
        end
        println(io)
        println(io, referencetext(p), ": ",
                allreferenced(p) ? "available in every spectrum" :
                "missing in some spectra")
        p.aligned[] && println(io, "Spectra shown aligned to the reference")
        println(io)
        println(io, "Multiplets (centres as measured, ppm):")
        for (m, mp) in enumerate(multiplets(t))
            ls = orderedlines(t, m)
            isempty(ls) && continue
            cs = filter(!isnothing, [linecentre(t, ls, i) for i in 1:n])
            span = isempty(cs) ? "no centre" : "$(fmt(minimum(cs))) to $(fmt(maximum(cs)))"
            println(io,
                    "  $(rpad(mp.label, 16)) $(length(ls)) line$(length(ls) == 1 ? "" : "s")" *
                    "   $(length(cs)) of $n spectra   $span")
        end
        isnothing(p.call) && return
        if isnothing(p.processing)
            writereproduce(io, folder, p.call; peaklist=joinpath(folder, "peaklist.csv"))
        else
            writereproduce(io, folder, p.call; lb=p.processing.lb[],
                           phase=joinpath(folder, "processing.csv"),
                           peaklist=joinpath(folder, "peaklist.csv"))
        end
        isfile(joinpath(folder, RECORD)) && p.reference isa IndirectReference &&
            println(io, "Pass reference= the 1H analysis, to reference indirectly again.")
        return
    end
    return path
end

"""
Copy the figures saved with the `f` key (`figure-*.pdf`) from the earlier folder that
`backupfolder` has just moved aside into the new one, so that they stay with the analysis.
"""
function copyfigures(dir)
    previous = rstrip(abspath(dir), ['/', '\\']) * "_previous"
    isdir(previous) || return dir
    for name in readdir(previous)
        occursin(r"^figure-.*\.pdf$", name) &&
            cp(joinpath(previous, name), joinpath(dir, name))
    end
    return dir
end

"Heading `heights.csv`, since a height is easily mistaken for an integral or an amplitude."
const HEIGHTNOTE = "Heights: the intensity at each line's maximum, divided by scale(spec) " *
                   "(ns x rg, and the concentration where recorded); heights, not " *
                   "integrals or fitted amplitudes"

"""
    saveanalysis(p, folder) -> String

Write the output folder: `summary.txt`, `peaklist.csv` (what was picked, which **Load**
reads back), `series.csv` (every line in every plane), the wide files `shifts.csv`,
`referencedshifts.csv` (once the reference is picked anywhere) and `heights.csv` (see
[`linetable`](@ref)), `processing.csv` for time-domain data, `stack.pdf` (the region shown
in the window), and in `peaks/` one CSV and one plot per multiplet. An existing folder is
moved aside first, and the figures saved into it with the `f` key are copied across.
Returns the folder.
"""
function saveanalysis(p::PeakTracker, folder::AbstractString)
    dir = backupfolder(folder)
    copyfigures(dir)
    comments = analysisinfo(p)
    writetable(joinpath(dir, "peaklist.csv"), comments, peaklisttable(p)...)
    header, rows = seriestable(p)
    writetable(joinpath(dir, "series.csv"), comments, header, rows)
    writetable(joinpath(dir, "shifts.csv"), [comments; "Shifts as measured"],
               linetable(p, q -> q.shift, "ppm")...)
    if isreferenced(p)
        writetable(joinpath(dir, "referencedshifts.csv"), [comments; referencetext(p)],
                   linetable(p, q -> referencedshift(p, q), "ppm")...)
    end
    writetable(joinpath(dir, "heights.csv"), [comments; HEIGHTNOTE],
               linetable(p, q -> q.height)...)
    isnothing(p.processing) ||
        writetable(joinpath(dir, "processing.csv"), comments,
                   processingtable(p.processing)...)
    save(joinpath(dir, "stack.pdf"), stackplot(p); backend=CairoMakie, pt_per_unit=1)

    t = tracks(p)
    for (m, mp) in enumerate(multiplets(t))
        isempty(orderedlines(t, m)) && continue
        base = joinpath(dir, "peaks", safename(mp.label))
        writetable(base * ".csv", comments, header, filter(row -> row[2] == mp.label, rows))
        save(base * ".pdf", multipletplot(p, m); backend=CairoMakie, pt_per_unit=1)
    end
    writesummary(joinpath(dir, "summary.txt"), p, dir)
    @info "Saved results to $dir"
    return dir
end
