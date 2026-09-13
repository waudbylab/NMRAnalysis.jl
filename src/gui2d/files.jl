# Peak list / results files (see writeresults! and readpeaklist!):
#   - lines beginning with '#' are comments (experiment metadata)
#   - an ordinary header row names the columns
#   - on read, only the label, x and y columns are used
# Hand-made lists may instead be a bare, header-less `label x y` per line.

function loadpeaks!(expt)
    file = pick_file(; filterlist="csv;list;peaks;txt;bak")
    file == "" && return

    @info "Loading peak file $file"
    isfitting = expt.isfitting[]
    if isfitting
        expt.isfitting[] = false
    end
    deleteallpeaks!(expt)
    readpeaklist!(expt, file)
    return expt.isfitting[] = isfitting
end

function saveresults!(expt)
    folder = joinpath(pwd(), expt.state[][:outputdir][])

    @info "Saving results to $folder"
    @async begin
        expt.state[][:mode][] = :fitting
        sleep(0.1) # allow time for mode change to be processed
    end
    @async begin # do saving in a separate task
        sleep(0.2) # allow time for mode change to be processed
        try
            # An existing folder is moved aside to <name>_previous rather than written
            # into, so a peak deleted since the last save leaves nothing behind - which is
            # what the stale-PDF sweep here used to be for.
            backupfolder(folder)
            writeresults!(expt, folder)
            writesummary(joinpath(folder, "summary.txt"), expt)
            save_peak_plots!(expt, folder)
            save_cluster_plots!(expt, folder)
            save_summary_plot!(expt, folder)
        catch e
            @error "Error saving results to $folder" exception = (e, catch_backtrace())
        finally
            # Always restore normal mode - otherwise a failure mid-save leaves the GUI stuck
            # on the salmon "fitting" background with no way to recover.
            GLMakie.activate!()
            expt.state[][:mode][] = :normal
        end
    end
end

"""
    splitfields(line)

Split a data line into fields, accepting either comma-separated (the format
written by [`writeresults!`](@ref)) or whitespace-separated values. Surrounding
whitespace on each field is stripped.
"""
function splitfields(line)
    fields = occursin(',', line) ? split(line, ',') : split(line)
    return strip.(fields)
end

"""
    stripunit(name) -> String

A column name with its parenthesised unit removed: `"x (ppm)"` → `"x"`. Column headers
carry units (see `docs/src/advanced/conventions.md`), so anything locating a column by name
strips them first.
"""
stripunit(name) = strip(replace(String(name), r"\s*\(.*\)\s*$" => ""))

"""
    headercolumns(filepath) -> Dict{String,Int} or nothing

The column positions of `filepath`'s header row, lowercased and with units stripped, or
`nothing` where the first non-comment line is data rather than a header (a hand-made
`label x y` list).
"""
function headercolumns(filepath::AbstractString)
    for line in eachline(filepath)
        sline = strip(line)
        (isempty(sline) || startswith(sline, '#')) && continue
        names = lowercase.(stripunit.(splitfields(sline)))
        return "label" in names ? Dict(n => i for (i, n) in enumerate(names)) : nothing
    end
    return nothing
end

"""
    parse_radius_comment!(expt, line)

If `line` records an X/Y fitting radius (as written by [`writeresults!`](@ref), e.g.
`# X radius / ppm: 0.04`), apply it to the experiment. Generic to all 2D experiments.
"""
function parse_radius_comment!(expt, line)
    mx = match(r"x\s*radius[^0-9]*([0-9]*\.?[0-9]+)"i, line)
    isnothing(mx) || setradius!(expt, :x, parse(Float64, mx.captures[1]))
    my = match(r"y\s*radius[^0-9]*([0-9]*\.?[0-9]+)"i, line)
    isnothing(my) || setradius!(expt, :y, parse(Float64, my.captures[1]))
    return
end

# Set a fitting radius. When the GUI is up, drive the slider (which keeps the display in sync
# via connect!); otherwise set the experiment observable directly.
function setradius!(expt, dim, value)
    g = get(expt.state[], :gui, nothing)
    g = g isa Observable ? g[] : g
    if g isa AbstractDict && haskey(g, :sgradii)
        slider = dim === :x ? g[:sgradii].sliders[1] : g[:sgradii].sliders[2]
        set_close_to!(slider, value)
    else
        (dim === :x ? expt.xradius : expt.yradius)[] = value
    end
    return
end

"""
    readpeaklist!(expt, filepath::AbstractString)

Read a peak list and add the peaks to `expt`. Only the `label`, `x` and `y`
columns are used — every other column (`resnum`, `resname`, `atom`, linewidths,
amplitudes, derived parameters …) is ignored, since the residue number and atom
are re-derived from the label.

Two layouts are accepted:
- **With a header row** (e.g. the program's own `results.csv`): the first
  non-comment line contains a `label` column name, and `x`/`y` are located by
  name, so column order and extra columns don't matter.
- **Without a header** (a hand-made list): values are read positionally as
  `label x y` from the first three columns, so no exact column names are needed.

Fields may be comma- or whitespace-separated. Lines beginning with `#` are
comments. A malformed line is skipped with a warning rather than aborting the
load — any labelling convention is tolerated.
"""
function readpeaklist!(expt, filepath::AbstractString)
    issparkylist(filepath) && return readsparkylist!(expt, filepath)

    columns = headercolumns(filepath)
    if !isnothing(columns)
        # A peak list written by this program carries a `plane` column, which is what
        # distinguishes one position per peak from a whole hand-tracked trajectory.
        haskey(columns, "plane") && return readtrackedpeaks!(expt, filepath, columns)
        # results.csv carries no positions: they are input, not output, and the fitted
        # ones are per plane. So look for the peak list beside it, then the series data.
        if !haskey(columns, "x") && !haskey(columns, "x[1]")
            for name in ("peaklist.csv", "series.csv")
                beside = joinpath(dirname(filepath), name)
                isfile(beside) && return readpeaklist!(expt, beside)
            end
            throw(ArgumentError("$filepath has no peak positions, and there is no " *
                                "peaklist.csv or series.csv beside it"))
        end
    end

    peak_count = 0
    colmap = nothing  # name => index, or nothing until established

    open(filepath) do f
        for (line_number, line) in enumerate(eachline(f))
            sline = strip(line)
            isempty(sline) && continue
            if startswith(sline, '#')
                parse_radius_comment!(expt, sline)  # restore fitting radii if recorded
                continue
            end

            fields = splitfields(sline)

            # Establish the column layout from the first non-comment line
            if isnothing(colmap)
                lower = lowercase.(stripunit.(fields))
                if "label" in lower
                    colmap = Dict(name => i for (i, name) in enumerate(lower))
                    continue  # header line, not data
                else
                    colmap = Dict("label" => 1, "x" => 2, "y" => 3)  # positional
                end
            end

            try
                length(fields) < 3 && error("insufficient fields")
                label = string(fields[colmap["label"]])

                # Moving-peak results store positions (and linewidths) per plane as
                # x[1], x[2], ...; restore them so a saved trajectory reloads intact. A
                # plain `label x y` list (or a fixed-peak results.csv) takes the single x/y.
                if !hasfixedpositions(expt) && haskey(colmap, "x[1]")
                    n = nslices(expt)
                    getcol(name, i) = parse(Float64, fields[colmap["$(name)[$i]"]])
                    xs = [getcol("x", i) for i in 1:n]
                    ys = [getcol("y", i) for i in 1:n]
                    addpeak!(expt, Point2f(xs[1], ys[1]), label)
                    peak = expt.peaks[][end]
                    setperplane!(peak, :x, xs)
                    setperplane!(peak, :y, ys)
                    if haskey(colmap, "r2x[1]")
                        setperplane!(peak, :R2x, [getcol("r2x", i) for i in 1:n])
                        setperplane!(peak, :R2y, [getcol("r2y", i) for i in 1:n])
                    end
                else
                    x = parse(Float64, fields[colmap["x"]])
                    y = parse(Float64, fields[colmap["y"]])
                    addpeak!(expt, Point2f(x, y), label)
                end
                peak_count += 1
            catch e
                @warn "Skipping line $line_number: $(e isa ErrorException ? e.msg : e)"
            end
        end
    end

    @debug "Added $peak_count peaks from $filepath"
    return peak_count
end

"""
    readtrackedpeaks!(expt, filepath, columns) -> Int

Restore peaks from a file holding one row per peak per plane - `peaklist.csv`, or a
`series.csv` read as a fallback. Each peak is added at its first row's position and the
whole trajectory is then set on it, so hand-tracking survives a save and reload.

A `peaklist.csv` row whose `plane` is blank means one position for every plane, so such a
peak is added with a single position rather than a trajectory; that is also what a
hand-made list or an imported Sparky list produces.
"""
function readtrackedpeaks!(expt, filepath::AbstractString, columns=headercolumns(filepath))
    (isnothing(columns) || !haskey(columns, "x") || !haskey(columns, "y")) &&
        throw(ArgumentError("$filepath has no label/x/y columns"))

    labels = String[]
    planes = Dict{String,Vector{Union{Nothing,Int}}}()
    positions = Dict{String,Vector{NTuple{4,Union{Nothing,Float64}}}}()
    for line in eachline(filepath)
        sline = strip(line)
        isempty(sline) && continue
        if startswith(sline, '#')
            parse_radius_comment!(expt, sline)
            continue
        end
        fields = splitfields(sline)
        lowercase(stripunit(fields[columns["label"]])) == "label" && continue  # header
        label = String(fields[columns["label"]])
        cell(name) = haskey(columns, name) && columns[name] ≤ length(fields) ?
                     tryparse(Float64, fields[columns[name]]) : nothing
        label in labels || push!(labels, label)
        push!(get!(positions, label, NTuple{4,Union{Nothing,Float64}}[]),
              (cell("x"), cell("y"), cell("r2x"), cell("r2y")))
        push!(get!(planes, label, Union{Nothing,Int}[]),
              haskey(columns, "plane") && columns["plane"] ≤ length(fields) ?
              tryparse(Int, fields[columns["plane"]]) : nothing)
    end

    # A row's plane is the plane it belongs to, not the position it happens to occupy in
    # the file: a list sorted by label or edited by hand can carry them out of order.
    for label in labels
        p = planes[label]
        all(!isnothing, p) && length(unique(p)) == length(p) || continue
        positions[label] = positions[label][sortperm(Int.(p))]
    end

    count = 0
    n = nslices(expt)
    for label in labels
        points = positions[label]
        xs = [p[1] for p in points]
        ys = [p[2] for p in points]
        (isempty(xs) || any(isnothing, xs) || any(isnothing, ys)) && continue
        addpeak!(expt, Point2f(xs[1], ys[1]), label)
        peak = expt.peaks[][end]
        # One row means one position for every plane, which `addpeak!` has already set.
        # A full trajectory is only restored when it matches this experiment's plane count;
        # a list carried over from a series of a different length seeds the first position
        # instead of silently mis-assigning the rest.
        if length(points) == n > 1
            setperplane!(peak, :x, Float64.(xs))
            setperplane!(peak, :y, Float64.(ys))
            r2x = [p[3] for p in points]
            r2y = [p[4] for p in points]
            if !any(isnothing, r2x) && !any(isnothing, r2y)
                setperplane!(peak, :R2x, Float64.(r2x))
                setperplane!(peak, :R2y, Float64.(r2y))
            end
        elseif length(points) > 1
            @warn "$label has $(length(points)) positions but this experiment has $n " *
                  "planes - using the first"
        end
        count += 1
    end
    @debug "Added $count peaks from $filepath"
    return count
end

"Sort peaks by residue number (positive ascending first, then unassigned)."
function sortedpeaks(expt)
    return sort(collect(expt.peaks[]);
                by=peak -> begin
                    r = extract_residue_number(peak.label[])
                    (r ≤ 0, abs(r))
                end)
end

"Format a post-fit parameter value/uncertainty, returning \"NA\" if absent."
function format_post(peak, key, which)
    haskey(peak.postparameters, key) || return "NA"
    val = getproperty(peak.postparameters[key], which)[][1]
    return string(to_value(val))
end

"""
    setperplane!(peak, param, values)

Write a per-plane parameter from a loaded results file, setting both the fitted value and the
initial value in every plane (so the loaded positions display immediately and seed any refit).
"""
function setperplane!(peak, param, values)
    haskey(peak.parameters, param) || return
    p = peak.parameters[param]
    for i in eachindex(values)
        p.value[][i] = values[i]
        p.initialvalue[][i] = values[i]
    end
    return
end

"""
    backup_file(filepath::AbstractString)

Create a backup of an existing file by appending '.old'.
"""
function backup_file(filepath::AbstractString)
    isfile(filepath) || return
    backup_path = filepath * ".old"
    @debug "Backing up $filepath to $backup_path"
    return mv(filepath, backup_path; force=true)
end

"""
    format_param(peak, param, slice, value_type) -> String

Format a parameter value, returning "NA" if not found.
"""
function format_param(peak, param, slice, value_type)
    haskey(peak.parameters, param) || return "NA"
    val = getproperty(peak.parameters[param], value_type)[][slice]
    return string(to_value(val)) # convert from Observables to plain values
end

"""
    save_cluster_plots!(expt, folder)

Save one zoomed contour plot per cluster of overlapping peaks (first plane only).
Axis limits are set to include only the peaks in that cluster plus padding.

Files are named `cluster_LABEL.pdf` (single peak) or `cluster_LABEL1-LABEL2.pdf`
(overlapping peaks).
"""
function save_cluster_plots!(expt, folder)
    peaks = expt.peaks[]
    clusters = expt.clusters[]
    isempty(clusters) && return

    CairoMakie.activate!()

    state = expt.state[]
    contourlevels = state[:gui][][:contourlevels]
    xlabel_str = "$(label(expt.specdata.nmrdata[1],F1Dim)) / ppm"
    ylabel_str = "$(label(expt.specdata.nmrdata[1],F2Dim)) / ppm"

    for cluster_idxs in clusters
        cluster = [peaks[i] for i in cluster_idxs]
        isempty(cluster) && continue

        # Bounding box from initial positions + radius padding
        padding = 1.5
        all_x = Float64[]
        all_y = Float64[]
        for peak in cluster
            pos = initialposition(peak)[]
            pts = pos isa AbstractVector ? pos : [pos]
            for p in pts
                push!(all_x, p[1])
                push!(all_y, p[2])
            end
        end
        max_xr = maximum(peak.xradius[] for peak in cluster)
        max_yr = maximum(peak.yradius[] for peak in cluster)
        lims = ((minimum(all_x) - padding * max_xr, maximum(all_x) + padding * max_xr),
                (minimum(all_y) - padding * max_yr, maximum(all_y) + padding * max_yr))

        # Build a filename from peak labels (truncate if very long)
        raw = join([peak.label[] for peak in cluster], "-")
        safe = replace(raw, r"[^\w\-]" => "_")
        safe = length(safe) > 64 ? safe[1:64] : safe

        # Always show only the first plane — avoids unwieldy grids for many-slice
        # experiments (CEST, relaxation series, etc.)
        fig = Figure(; size=(350, 350))
        ax = Axis(fig[1, 1];
                  xlabel=xlabel_str, ylabel=ylabel_str,
                  xreversed=true, yreversed=true, limits=lims)
        heatmap!(ax, expt.specdata.x[1], expt.specdata.y[1],
                 expt.specdata.mask[][1];
                 colormap=[:white, :lightgoldenrod1], colorrange=(0, 1))
        contour!(ax, expt.specdata.x[1], expt.specdata.y[1],
                 expt.specdata.z[1];
                 levels=contourlevels, color=bicolours(:grey50, :lightblue))
        contour!(ax, expt.specdata.x[1], expt.specdata.y[1],
                 expt.specdata.zfit[][1];
                 levels=contourlevels, color=bicolours(:orangered, :dodgerblue))
        for peak in cluster
            pt = Point2f(peak.parameters[:x].value[][1],
                         peak.parameters[:y].value[][1])
            scatter!(ax, [pt]; markersize=10, marker=:x, color=:black)
            text!(ax, [pt]; text=[peak.label[]], fontsize=12,
                  offset=(6, 0), align=(:left, :center), color=:black)
        end

        save(joinpath(folder, "cluster_$(safe).pdf"), fig)
    end

    return GLMakie.activate!()
end

"""
    save_summary_plot!(expt, folder)

Write `summary.pdf` to `folder` using the experiment's default summary parameter.
Skipped when there are no peaks.
"""
function save_summary_plot!(expt, folder)
    isempty(expt.peaks[]) && return
    CairoMakie.activate!()
    try
        fig = summaryplot(expt)
        save(joinpath(folder, "summary.pdf"), fig)
    finally
        GLMakie.activate!()
    end
end
# ---- Sparky peak lists --------------------------------------------------------
# Sparky writes a whitespace-delimited list with an `Assignment w1 w2` header and, often,
# further columns (data height, volume, notes) that are ignored here. It carries one
# position per peak and nothing else: no trajectory, no fitting radii, no uncertainties.
# That is why it is an import format rather than the peak list this program writes - see
# `peaklisttable` in output.jl.

"""
    issparkylist(filepath) -> Bool

Whether `filepath` looks like a Sparky peak list: a whitespace-delimited file whose first
non-comment line names an `Assignment` column.
"""
function issparkylist(filepath::AbstractString)
    for line in eachline(filepath)
        sline = strip(line)
        (isempty(sline) || startswith(sline, '#')) && continue
        occursin(',', sline) && return false
        return lowercase(first(split(sline))) == "assignment"
    end
    return false
end

"""
    readsparkylist!(expt, filepath) -> Int

Import a Sparky peak list, adding its peaks to `expt`.

Sparky names its dimensions `w1`, `w2` in the spectrum's own order, which for a ¹⁵N-HSQC
conventionally puts the indirect dimension first - the opposite way round from this
program, whose `x` is the direct dimension. Reading the columns in order would therefore
transpose every peak, so which column goes to which axis is decided by
[`sparkyaxisorder`](@ref) from where the shifts actually fall, and only falls back to the
convention when that cannot tell.

The assignment string becomes the peak label unchanged: Sparky's `G10N-G10H` and a plain
`G10` both survive a round trip, and `parse_label` already derives the residue number and
atom from whatever form it is given.
"""
function readsparkylist!(expt, filepath::AbstractString)
    entries = Tuple{String,Float64,Float64}[]
    for line in eachline(filepath)
        sline = strip(line)
        (isempty(sline) || startswith(sline, '#')) && continue
        fields = split(sline)
        length(fields) < 3 && continue
        lowercase(fields[1]) == "assignment" && continue        # header
        w1 = tryparse(Float64, fields[2])
        w2 = tryparse(Float64, fields[3])
        (isnothing(w1) || isnothing(w2)) && continue
        push!(entries, (String(fields[1]), w1, w2))
    end
    isempty(entries) && throw(ArgumentError("$filepath contains no Sparky peaks"))

    swapped = sparkyaxisorder(expt, entries)
    for (label, w1, w2) in entries
        x, y = swapped ? (w1, w2) : (w2, w1)
        addpeak!(expt, Point2f(x, y), label)
    end
    @debug "Added $(length(entries)) peaks from Sparky list $filepath"
    return length(entries)
end

"""
    sparkyaxisorder(expt, entries) -> Bool

Whether a Sparky list's `w1` belongs on this experiment's `x` (direct) axis rather than its
`y` axis - `true` for a list written the opposite way round from the usual convention.

Decided by counting how many peaks land inside both axes' actual chemical-shift ranges each
way round, because that is the thing that is really being asked and it needs no metadata
beyond the spectrum itself. A tie (including the degenerate case of two dimensions covering
similar ranges, as in a NOESY) keeps Sparky's convention of `w1` on the indirect axis.
"""
function sparkyaxisorder(expt, entries)
    spec = expt.specdata.nmrdata[1]
    xrange = extrema(data(spec, F1Dim))
    yrange = extrema(data(spec, F2Dim))
    inside(v, range) = range[1] ≤ v ≤ range[2]
    conventional = count(e -> inside(e[3], xrange) && inside(e[2], yrange), entries)
    swapped = count(e -> inside(e[2], xrange) && inside(e[3], yrange), entries)
    swapped > conventional &&
        @info "Sparky list appears transposed (w1 on the direct axis) - reading it that way"
    return swapped > conventional
end
