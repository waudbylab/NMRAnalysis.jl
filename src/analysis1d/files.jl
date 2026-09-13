# Output files, written to the layout and column rules in
# `docs/src/advanced/conventions.md`:
#
#   out/
#     summary.txt        (written by `saveresults` in gui.jl)
#     regionlist.csv     what the user picked; what Load reads back
#     results.csv        one row per region: everything reported for it
#     series.csv         the measurements, one row per region per plane
#     fit.pdf
#     regions/<label>.csv, regions/<label>.pdf

# ---- value formatting ---------------------------------------------------------
# `csvcolumn`, `csvcolumns`, `csvvalue`, `safename`, `backupfile` and the underlying
# `writetable` are shared with the other modules - see `src/output.jl`.

"""Value and uncertainty cells for a parameter that may carry neither."""
function valueerr(params, key::Symbol)
    haskey(params, key) || return ("NA", "NA")
    v = params[key]
    v isa Measurement || return (csvvalue(v), "NA")
    return (csvvalue(Measurements.value(v)), csvvalue(Measurements.uncertainty(v)))
end

# ---- provenance ---------------------------------------------------------------

"""
    experimentinfo(expt, [dataset]) -> String

Multi-line description of the experiment, written as the comment header of every CSV.
Generic across every experiment - the analysis name, where the data came from, how many
spectra and which variables are arrayed are all already known to the framework - so
experiments need not override it (GUI2D's equivalent is written out per experiment).

The dataset is passed explicitly because the one the experiment was *constructed* with
still carries the noise position it started at, where the GUI's live dataset carries
wherever the user dragged the marker to. Provenance has to describe the data as analysed.
"""
function experimentinfo(expt::Experiment1D, ds::Dataset1D=dataset(expt))
    io = IOBuffer()
    println(io, "NMRAnalysis.jl $(pkgversion(@__MODULE__))")
    println(io, "Analysis: $(windowtitle(expt))")
    for src in unique(sources(ds))
        isempty(src) || println(io, "Source: $src")
    end
    println(io, "Number of spectra: $(nplanes(ds))")
    if !isempty(ds.planes.vars)
        println(io, "Arrayed variables: $(join(keys(first(ds.planes.vars)), ", "))")
    end
    println(io, "Noise position / ppm: $(round(ds.noisecenter; digits=4))")
    return String(take!(io))
end

# ---- results.csv --------------------------------------------------------------

"""
    resultkeys(expt) -> Vector{Symbol}

The grouping coordinates that, with the region label, identify a row of `results.csv`:
`(:which,)` for TRACT, `(:run,)` for a multi-run kinetics series, none for a single-series
experiment.
"""
resultkeys(expt::Experiment1D) = collect(Symbol, groupcols(expt))

"""
    parameternames(results) -> Vector{Symbol}

Every parameter name present across `results`, in order of first appearance. The union, so
the table stays rectangular when one region derived something another did not.
"""
function parameternames(results)
    names = Symbol[]
    for r in results, name in keys(r.parameters)
        name in names || push!(names, name)
    end
    return names
end

"""
    resultstable(expt, results) -> (header, rows)

Column names and rows for `results.csv`: one row per region, carrying every parameter
reported for it, with the experiment's [`primaryparam`](@ref) first. A region measured
under several conditions keeps them all on that row under the names [`seriesname`](@ref)
gave them, so TRACT writes `R_trosy` and `R_anti` beside the `tauc` computed from the two.

The region bounds are not here: they are what the user picked, and live in
`regionlist.csv` (see [`regionlisttable`](@ref)).
"""
function resultstable(expt::Experiment1D, results)
    names = parameternames(results)
    primary = primaryparam(expt)
    i = findfirst(==(primary), names)
    isnothing(i) || pushfirst!(names, popat!(names, i))

    header = ["label"]
    for name in names
        # the unit belongs to the quantity, so `R_trosy` is looked up as `R`
        append!(header, collect(csvcolumns(name, paramunit(expt, baseparam(name)))))
    end

    rows = Vector{String}[]
    for r in results
        row = [r.region]
        for name in names
            append!(row, collect(valueerr(r.parameters, name)))
        end
        push!(rows, row)
    end
    return header, rows
end

# ---- series.csv ---------------------------------------------------------------

"""
    fitvalues(series) -> Vector{Float64}

The model evaluated at the measured coordinates, for the `I_fit` column - not on a fine
grid, so that a residual is a subtraction. `NaN` (written `NA`) where nothing was fitted.
"""
fitvalues(s::SeriesResult) = fitvalues(s.model, s)
fitvalues(::SeriesModel, s::SeriesResult) = fill(NaN, length(s.x))
function fitvalues(m::CurveFitModel, s::SeriesResult)
    isempty(s.coefficients) && return fill(NaN, length(s.x))
    return m.func(s.x, Measurements.value.(s.coefficients))
end

"""
    seriestable(expt, dataset, results) -> (header, rows)

Column names and rows for `series.csv`: the measurements themselves, long, one row per
region per plane.

`source` and `plane` are always written, even when constant. They are provenance, so
whatever physical variable distinguishes several datasets (TRACT's `which`, a
concentration) gets its own coordinate column as well. `plane` is load-bearing where the
planes share one file: every row of a pseudo-2D experiment has the same `source`.
"""
function seriestable(expt::Experiment1D, ds::Dataset1D, results)
    keycols = resultkeys(expt)
    axis = fitaxis(expt)
    header = ["source", "label", "plane"]
    append!(header, [csvcolumn(k, coordinateunit(expt, k)) for k in keycols])
    push!(header, csvcolumn(axis, coordinateunit(expt, axis)))
    append!(header, ["I", "I_err", "I_fit"])

    src = sources(ds)
    rows = Vector{String}[]
    for r in results, s in r.series
        yfit = fitvalues(s)
        for i in eachindex(s.x)
            plane = i ≤ length(s.planes) ? s.planes[i] : 0
            row = [1 ≤ plane ≤ length(src) ? src[plane] : ds.label, r.region,
                   csvvalue(plane)]
            append!(row, [csvvalue(get(s.group, k, nothing)) for k in keycols])
            push!(row, csvvalue(s.x[i]))
            append!(row,
                    [csvvalue(Measurements.value(s.y[i])),
                     csvvalue(Measurements.uncertainty(s.y[i])),
                     csvvalue(yfit[i])])
            push!(rows, row)
        end
    end
    return header, rows
end

# ---- writing ------------------------------------------------------------------

"""
    writeresults!(expt, dataset, results, regions, folder) -> String

Write the whole result set into `folder`: `regionlist.csv` (what was picked), `series.csv`,
`results.csv` where anything is reported per region, and one `regions/<label>.csv` per
region holding that region's own rows of `series.csv`. Returns the path of the series
file.
"""
function writeresults!(expt::Experiment1D, ds::Dataset1D, results, regs,
                       folder::AbstractString)
    comments = split(experimentinfo(expt, ds), '\n')

    # What the user picked, kept apart from what the fit produced, so a region list can be
    # reused on another dataset (and is what `readregions!` reads back).
    defaultwidth = defaultregionwidth(first(ds.planes.traces).δ)
    writetable(joinpath(folder, "regionlist.csv"), comments,
               regionlisttable(regs, ds.noisecenter, defaultwidth)...)

    # An experiment that reports nothing per region (kinetics) would get a results.csv of
    # labels and no values, so it gets none - see docs/src/advanced/conventions.md.
    isempty(parameternames(results)) ||
        writetable(joinpath(folder, "results.csv"), comments,
                   resultstable(expt, results)...)

    header, rows = seriestable(expt, ds, results)
    filepath = writetable(joinpath(folder, "series.csv"), comments, header, rows)

    # per-region files, beside the per-region plots `saveresults` writes
    labelcol = findfirst(==("label"), header)
    for reg in regs
        writetable(joinpath(folder, "regions", "$(safename(reg.label)).csv"), comments,
                   header, filter(row -> row[labelcol] == reg.label, rows))
    end

    return filepath
end

# ---- regionlist.csv: the user's input -----------------------------------------

"""
    NOISE_LABEL

The label of the row in `regionlist.csv` that records the noise position rather than a
signal region. Reserved: a signal region of this name would be read back as the noise
marker.
"""
const NOISE_LABEL = "noise"

"""
    regionlisttable(regions, noisecentre, defaultwidth) -> (header, rows)

Column names and rows for `regionlist.csv`: what the user picked, as opposed to what the
fit produced. One row per region, plus one for the noise position.

The noise marker has a position but no width of its own - the width used to estimate a
region's uncertainty always matches that region's own (see [`integrate`](@ref)) - so it
is written as a region `defaultwidth` wide, or as wide as the widest signal region where
that is wider. Reading its centre back is what matters; the width is there so the row means
something on its own and so the file needs no separate convention for it.
"""
function regionlisttable(regs, noisecentre, defaultwidth)
    header = ["label", csvcolumn("lo", "ppm"), csvcolumn("hi", "ppm")]
    rows = [[r.label, csvvalue(r.lo), csvvalue(r.hi)] for r in regs]
    w = maximum([defaultwidth; [width(r) for r in regs]])
    push!(rows, [NOISE_LABEL, csvvalue(noisecentre - w / 2), csvvalue(noisecentre + w / 2)])
    return header, rows
end

"""
    readregions!(state, filepath) -> Int

Restore a saved region list, replacing whatever is currently set, and return the number of
signal regions read. The row labelled `noise` sets the noise position from its centre;
every other row is a region.

Written for `regionlist.csv` but deliberately tolerant, so that a hand-made list works and
so does a `results.csv` from before the region bounds moved out of it: any file with a
label column and `lo`/`hi` columns is read, and a `# Noise position / ppm:` comment is
honoured where there is no `noise` row.
"""
function readregions!(state, filepath::AbstractString)
    isfile(filepath) || throw(ArgumentError("no such region list: $filepath"))
    regs = Region[]
    seen = Set{String}()
    colmap = nothing
    labelcol = nothing
    for line in eachline(filepath)
        sline = strip(line)
        isempty(sline) && continue
        if startswith(sline, '#')
            m = match(r"Noise position / ppm:\s*(\S+)", sline)
            isnothing(m) || (state[:noisec][] = parse(Float64, m.captures[1]))
            continue
        end
        fields = strip.(split(sline, ','))
        if isnothing(colmap)
            # Units are part of the header (`lo (ppm)`), so columns are located by their
            # name with any parenthesised unit stripped.
            colmap = Dict(lowercase(strip(replace(name, r"\s*\(.*\)$" => ""))) => i
                          for (i, name) in enumerate(fields))
            labelcol = get(colmap, "label", get(colmap, "region", nothing))
            (isnothing(labelcol) || !haskey(colmap, "lo") || !haskey(colmap, "hi")) &&
                throw(ArgumentError("$filepath has no label/lo/hi columns"))
            continue
        end
        label = String(fields[labelcol])
        (isempty(label) || label in seen) && continue
        lo, hi = fields[colmap["lo"]], fields[colmap["hi"]]
        (lo == "NA" || hi == "NA") && continue
        push!(seen, label)
        if lowercase(label) == NOISE_LABEL
            state[:noisec][] = (parse(Float64, lo) + parse(Float64, hi)) / 2
        else
            push!(regs, Region(label, parse(Float64, lo), parse(Float64, hi)))
        end
    end
    isempty(regs) && throw(ArgumentError("$filepath contains no regions"))
    state[:regions][] = regs
    state[:active][] = 1
    return length(regs)
end

# ---- the call that produced an analysis ---------------------------------------

"""
    AnalysisCall

The call that produced an analysis, recorded by the entry point so `summary.txt` can print
a line that repeats it. `reproducible` is false where one of the positional arguments has
no Julia literal - a spectrum passed as an `NMRData` rather than as a path - in which case
no such line is printed at all rather than a misleading one.
"""
struct AnalysisCall
    func::String
    args::Vector{String}
    kwargs::Vector{Pair{Symbol,String}}
    reproducible::Bool
end

"""
    analysiscall(func, args...; kwargs...) -> AnalysisCall

Record a call. Keyword arguments whose values have no literal (a `SeriesModel` object, say)
are left out, since the point of the record is that it can be pasted back into Julia.
"""
function analysiscall(func::AbstractString, args...; kwargs...)
    rendered = Pair{Symbol,String}[]
    for (name, value) in pairs(kwargs)
        text = callvalue(value)
        isnothing(text) || push!(rendered, name => text)
    end
    args_ = [callvalue(a) for a in args]
    return AnalysisCall(String(func), String[something(a, "spec") for a in args_],
                        rendered, !any(isnothing, args_))
end

"""
    callvalue(x) -> String or nothing

`x` as a Julia literal, or `nothing` where it has none. Vectors of numbers are written out
in full: a gradient ramp or a delay list is exactly what has to be repeated, so abbreviating
it would defeat the purpose.
"""
callvalue(x::Union{Real,Bool,Symbol,AbstractString}) = repr(x)
callvalue(x::AbstractVector{<:Real}) = repr(collect(x))
callvalue(x::NamedTuple) = repr(x)
callvalue(::Nothing) = nothing
callvalue(x) = nothing

"""
    callstring(call, regions, noisecentre) -> String or nothing

The `Reproduce:` line for `summary.txt`: the recorded call with the region actually used
appended as an `integration` triple, so pasting it back repeats the analysis without the
window opening.

Only a single-region analysis can be repeated in one call, the integration triple carrying
one region and one noise position. With several regions the call is still printed, without
them, and [`writesummary`](@ref) adds a line pointing at the **Load** button, which restores
every region *and* the noise position from `results.csv`.
"""
function callstring(call::AnalysisCall, regs, noisecentre)
    call.reproducible || return nothing
    kwargs = copy(call.kwargs)
    if length(regs) == 1
        r = only(regs)
        push!(kwargs,
              :integration =>
                  "(peakppm=$(round(centre(r); digits=3)), " *
                  "noiseppm=$(round(noisecentre; digits=3)), " *
                  "ppmwidth=$(round(width(r); digits=3)))")
    end
    io = IOBuffer()
    print(io, call.func, "(", join(call.args, ", "))
    if isempty(kwargs)
        print(io, ")")
    else
        isempty(call.args) || print(io, ";")
        println(io)
        pad = " "^(length(call.func) + 1)
        for (i, (name, value)) in enumerate(kwargs)
            print(io, pad, name, "=", value, i == length(kwargs) ? ")" : ",\n")
        end
    end
    return String(take!(io))
end

# ---- summary.txt --------------------------------------------------------------

"""
    writesummary(filepath, expt, dataset, results, regions, call=nothing)

Write `summary.txt`: the human-readable record of the analysis, and the one output where
numbers are rounded for reading rather than written at full precision.
"""
function writesummary(filepath::AbstractString, expt::Experiment1D, ds::Dataset1D, results,
                      regs, call=nothing)
    backupfile(filepath)
    open(filepath, "w") do f
        for line in split(experimentinfo(expt, ds), '\n')
            isempty(strip(line)) && continue
            println(f, line)
        end
        println(f)
        println(f, "Integration regions / ppm:")
        for r in regs
            println(f, "  $(r.label): $(round(r.lo; digits=4)) to $(round(r.hi; digits=4))")
        end
        println(f)
        for r in regs
            print(f, summarytext(expt, results, r.label))
        end
        isnothing(call) && return nothing
        text = callstring(call, regs, ds.noisecenter)
        isnothing(text) && return nothing
        println(f, "Reproduce:")
        for line in split(text, '\n')
            println(f, "    ", line)
        end
        if length(regs) != 1
            println(f)
            println(f,
                    "Then press Load in the analysis window to restore the " *
                    "$(length(regs)) regions and the noise position from results.csv.")
        end
        return nothing
    end
    return filepath
end
