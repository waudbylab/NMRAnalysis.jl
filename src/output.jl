# Shared output primitives.
#
# The column and value rules from `docs/src/advanced/conventions.md`, in one place, so that
# every analysis module writes the same kind of file. Deliberately only the *formatting* and
# *writing* is shared: building the tables needs to know what an experiment is, so that
# stays in each module (`analysis1d/files.jl`, `exchange1d/files.jl`, ...).

"""
    csvcolumn(name, unit) -> String

A CSV column header: the ASCII parameter or coordinate name, with its unit in parentheses
where it has one (`"R (s-1)"`, `"A"`). Units are ASCII so that a file survives being opened
on any machine, and a reader locates a column by stripping the parenthesised part.
"""
csvcolumn(name, unit) = isempty(unit) ? string(name) : "$(name) ($(unit))"

"""
    csvcolumns(name, unit) -> (String, String)

Value and uncertainty column headers for one parameter. The uncertainty repeats the unit,
so the pair is symmetrical for anything reading them.
"""
csvcolumns(name, unit) = (csvcolumn(name, unit), csvcolumn("$(name)_err", unit))

"""
    csvvalue(x) -> String

One cell, at full precision - these are machine files, and rounding is irreversible.
`"NA"` for anything absent or not applicable, which is distinct from a blank key (a blank
key means the row applies to every value of it).
"""
csvvalue(::Nothing) = "NA"
csvvalue(::Missing) = "NA"
csvvalue(x::Real) = isfinite(x) ? string(x) : "NA"
csvvalue(x) = string(x)

"""
    safename(label) -> String

A label made safe to use as a filename, so that an entity the user named anything at all
still writes to `<name>.csv` beside its plot without escaping its folder.
"""
function safename(label::AbstractString)
    name = replace(String(label), r"[^A-Za-z0-9._-]" => "_")
    return isempty(name) ? "unnamed" : name
end

"""
    sanitizelabel(label) -> String

Strip commas from a user-entered region or peak label. Labels are written unescaped into
CSV cells (see [`writetable`](@ref)), so a comma in one would corrupt the row.
"""
sanitizelabel(label::AbstractString) = replace(String(label), "," => "")

"""
    backupfile(filepath)

Rename an existing file to `<name>.bak`, so a save never silently destroys the last one.
"""
function backupfile(filepath::AbstractString)
    isfile(filepath) || return nothing
    return mv(filepath, filepath * ".bak"; force=true)
end

"""
    backupfolder(folder) -> String

Move an existing output folder aside to `<folder>_previous` and return the (now empty)
path, so that a save never silently destroys the last one and never leaves stale files from
it behind. An earlier `_previous` is replaced.

Folder-level rather than file-level, because the stale files are the problem: a region or a
peak deleted between one save and the next would otherwise leave its plot and its data in
place, looking like part of the current result.
"""
function backupfolder(folder::AbstractString)
    path = rstrip(abspath(folder), ['/', '\\'])
    cwd = rstrip(abspath(pwd()), ['/', '\\'])
    # An output-folder box left empty resolves to the working directory itself, and the
    # move below would take the running session's directory with it.
    (path == cwd || startswith(cwd, path * "/")) &&
        throw(ArgumentError("refusing to save into $path: it is, or contains, the working directory"))
    isdir(folder) && mv(folder, path * "_previous"; force=true)
    mkpath(folder)
    return folder
end

"""
    writetable(filepath, comments, header, rows) -> String

Write one CSV: each line of `comments` as a `#` comment (the experiment description), then
the `header` row, then `rows`. Any existing file is backed up first. Returns the path.
"""
function writetable(filepath::AbstractString, comments, header, rows)
    backupfile(filepath)
    mkpath(dirname(filepath))
    open(filepath, "w") do f
        for line in comments
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
