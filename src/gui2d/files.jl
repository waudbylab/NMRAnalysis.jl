# Peak list / results files (see writeresults! and readpeaklist!):
#   - lines beginning with '#' are comments (experiment metadata)
#   - an ordinary header row names the columns
#   - on read, only the label, x and y columns are used
# Hand-made lists may instead be a bare, header-less `label x y` per line.

function loadpeaks!(expt)
    file = pick_file(; filterlist="csv;peaks;txt;old")
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
    folder = pick_folder()
    folder == "" && return

    @info "Saving results to $folder"
    @async begin
        expt.state[][:mode][] = :fitting
        sleep(0.1) # allow time for mode change to be processed
    end
    @async begin # do fitting in a separate task
        sleep(0.2) # allow time for mode change to be processed

        # save all peak positions, linewidths, amplitudes and derived
        # parameters to a single results file
        writeresults!(expt, folder)

        # remove any files peak_XXX.pdf from output folder before saving new plots
        for file in readdir(folder)
            if occursin(r"^peak_.*\.pdf$", file)
                rm(joinpath(folder, file))
            end
        end
        save_peak_plots!(expt, folder)

        # save current contour plot figure
        save_contour_plot!(expt, folder)

        expt.state[][:mode][] = :normal
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
    peak_count = 0
    colmap = nothing  # name => index, or nothing until established

    open(filepath) do f
        for (line_number, line) in enumerate(eachline(f))
            sline = strip(line)
            isempty(sline) && continue
            startswith(sline, '#') && continue

            fields = splitfields(sline)

            # Establish the column layout from the first non-comment line
            if isnothing(colmap)
                lower = lowercase.(fields)
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
                x = parse(Float64, fields[colmap["x"]])
                y = parse(Float64, fields[colmap["y"]])

                addpeak!(expt, Point2f(x, y), label)
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
    writeresults!(expt, folder) -> String

Write all peak results to a single `results.csv` in `folder`. Each row is one
peak with identity (`label`, `resnum`, `resname`, `atom`), positions (`x`, `y`),
linewidths (`R2x`, `R2y`), per-plane amplitudes (`amp[1]`, `amp[2]`, …) and any
derived parameters, each value immediately followed by its `_err` uncertainty
column. Experiment metadata is written as `#`-comment lines above an ordinary
(uncommented) header row, so the file opens directly in spreadsheets and via
`pandas.read_csv(comment="#")`.
"""
function writeresults!(expt, folder)
    filepath = joinpath(folder, "results.csv")
    backup_file(filepath)

    header, rows = resultstable(expt)
    open(filepath, "w") do f
        for line in split(experimentinfo(expt), '\n')
            isempty(strip(line)) && continue
            println(f, "# ", line)
        end
        println(f, join(header, ","))
        for row in rows
            println(f, join(row, ","))
        end
    end
    return filepath
end

"Sort peaks by residue number (positive ascending first, then unassigned)."
function sortedpeaks(expt)
    return sort(collect(expt.peaks[]); by=peak -> begin
                    r = extract_residue_number(peak.label[])
                    (r ≤ 0, abs(r))
                end)
end

"""
    resultstable(expt) -> (header, rows)

Build the column-name `header` and the `rows` (each a vector of strings) for the
results file. Derived (post-fit) parameters are appended with the experiment's
primary parameter first (see `primaryparam`).
"""
function resultstable(expt)
    n = nslices(expt)
    peaks = sortedpeaks(expt)

    # derived parameter keys, primary result first
    derivedkeys = Symbol[]
    if !isempty(peaks)
        allkeys = collect(keys(first(peaks).postparameters))
        prim = primaryparam(expt)
        derivedkeys = prim in allkeys ? [prim; filter(!=(prim), allkeys)] : allkeys
    end

    header = ["label", "resnum", "resname", "atom",
              "x", "x_err", "y", "y_err", "R2x", "R2x_err", "R2y", "R2y_err"]
    for i in 1:n
        append!(header, ["amp[$i]", "amp[$i]_err"])
    end
    for k in derivedkeys
        append!(header, [string(k), "$(k)_err"])
    end

    rows = Vector{String}[]
    for peak in peaks
        lbl = parse_label(peak.label[])
        row = [peak.label[], string(lbl.resnum),
               lbl.onelettercode == '?' ? "" : string(lbl.onelettercode),
               lbl.atom]
        for p in (:x, :y, :R2x, :R2y)
            push!(row, format_param(peak, p, 1, :value))
            push!(row, format_param(peak, p, 1, :uncertainty))
        end
        for i in 1:n
            push!(row, format_param(peak, :amp, i, :value))
            push!(row, format_param(peak, :amp, i, :uncertainty))
        end
        for k in derivedkeys
            push!(row, format_post(peak, k, :value))
            push!(row, format_post(peak, k, :uncertainty))
        end
        push!(rows, row)
    end
    return header, rows
end

"Format a post-fit parameter value/uncertainty, returning \"NA\" if absent."
function format_post(peak, key, which)
    haskey(peak.postparameters, key) || return "NA"
    val = getproperty(peak.postparameters[key], which)[][1]
    return string(to_value(val))
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

function save_contour_plot!(expt, folder)
    CairoMakie.activate!()

    state = expt.state[]
    rect = state[:gui][][:axcontour].finallimits[]
    x0, y0 = rect.origin
    dx, dy = rect.widths
    lims = ((x0, x0 + dx), (y0, y0 + dy))

    for i in 1:nslices(expt)
        f = Figure()
        ax = Axis(f[1, 1];
                  xlabel="$(label(expt.specdata.nmrdata[1],F1Dim)) chemical shift (ppm)",
                  ylabel="$(label(expt.specdata.nmrdata[1],F2Dim)) chemical shift (ppm)",
                  xreversed=true, yreversed=true, limits=lims,
                  title=slicelabel(expt, i))
        heatmap!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.mask[][i];
                 colormap=[:white, :lightgoldenrod1],
                 colorrange=(0, 1))
        contour!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.z[i];
                 levels=state[:gui][][:contourlevels],
                 color=bicolours(:grey50, :lightblue))
        contour!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.zfit[][i];
                 levels=state[:gui][][:contourlevels],
                 color=bicolours(:orangered, :dodgerblue))
        scatter!(ax, state[:positions]; markersize=10, marker=:x, color=:black)
        text!(ax, state[:positions]; text=state[:labels],
              fontsize=14,
              # font=:bold,
              offset=(8, 0),
              align=(:left, :center),
              color=:black)

        save(joinpath(folder, "spectrum-fit-$i.pdf"), f)
    end

    return GLMakie.activate!()
end