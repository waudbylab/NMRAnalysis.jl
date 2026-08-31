last_folder = Ref{String}("")

"""
    select_expts(directory_path=""; experiment_type_filter="")

Interactively select NMR experiment folders from a given directory path.

# Arguments
- `directory_path::String`: Path to start folder selection from. If empty, uses last selected folder or current directory.
- `experiment_type_filter`: Optional filter to apply for experiment types.

# Returns
- `Vector{String}`: Array of selected experiment folder paths. Returns empty array if no selection is made.
"""
function select_expts(exptno::Integer; kwargs...)
    return select_expts(string(exptno); kwargs...)
end

function select_expts(directory_path=""; experiment_type_filter="")
    if isempty(directory_path)
        directory_path = last_folder[]
    end
    if isempty(directory_path)
        directory_path = pwd()
    end

    # If the given path is an experiment, return it in an array
    if isexpt(directory_path)
        return [directory_path]
    end

    if foldercontainsexperiments(directory_path)
        return pick_expt_folders(directory_path; experiment_type_filter)
    end

    # Use native file picker to select a folder from the given directory path
    folder = pick_folder(directory_path)

    # Return an empty array if no folder is selected
    if isempty(folder)
        return []
    end
    last_folder[] = folder

    # If the selected folder is an experiment, return it in an array
    if isexpt(folder)
        return [folder]
    else
        # Otherwise, pick experiment folders within the selected folder
        return pick_expt_folders(folder; experiment_type_filter)
    end
end

function foldercontainsexperiments(directory_path)
    numbered_dirs = filter(ispath, readdir(directory_path; join=true))
    for dir in numbered_dirs
        if isexpt(dir)
            return true
        end
    end
    return false
end

function types_and_features(filename)
    try
        expt = loadnmr(filename)
        hasannotations(expt) || return (String[], String[])
        expt_types = annotations(expt, :experiment_type)
        features = annotations(expt, :features)
        return (expt_types, features)
    catch
        return (String[], String[])
    end
end

"""
    pick_expt_folders(directory_path; experiment_type_filter="") -> Vector{String}
"""
function pick_expt_folders(directory_path; experiment_type_filter="")
    # Get all experiment directories and their titles
    numbered_dirs = filter(ispath, readdir(directory_path; join=true))
    filter!(isexpt, numbered_dirs)

    # Get types and features for all experiments (only load once per file)
    annotations_dict = Dict{String,Tuple{Vector{String},Vector{String}}}()
    for dir in numbered_dirs
        annotations_dict[dir] = types_and_features(dir)
    end

    # Experiment type filter function
    function type_filter_func(filename)
        experiment_type_filter == "" && return true
        types, _ = annotations_dict[filename]
        return experiment_type_filter in types
    end
    filter!(type_filter_func, numbered_dirs)

    # Build menu options
    options = String[]
    folder_paths = String[]
    for dir in numbered_dirs
        folder_num = basename(dir)

        title_path = joinpath(dir, "pdata", "1", "title")
        title_content = strip(readline(title_path))

        # Get experiment types and features from cache
        types, features = annotations_dict[dir]

        # Build the display string
        display_string = "$(folder_num): $(title_content)"

        # Add types and features if they exist
        if !isempty(types) || !isempty(features)
            annotations_str = "("
            if !isempty(types)
                annotations_str *= join(types, ", ")
            end
            if !isempty(features)
                if !isempty(types)
                    annotations_str *= "; "
                end
                annotations_str *= join(features, ", ")
            end
            annotations_str *= ")"
            display_string *= " " * annotations_str
        end

        push!(options, display_string)
        push!(folder_paths, dir)
    end

    # Handle case where no folders match filters
    if isempty(options)
        error("No r1rho experiments found in '$directory_path'.")
    end

    # Sort both arrays based on folder numbers
    perm = sortperm(options; by=x -> parse(Int, split(x, ':')[1]))
    options = options[perm]
    folder_paths = folder_paths[perm]

    # Create and display the multi-select menu
    menu = MultiSelectMenu(options)
    choices = request("Select experiments to analyse:", menu)

    # Convert Set to Vector for indexing and return the full paths
    return folder_paths[collect(choices)]
end

"""
    isexpt(directory) -> Bool

Check if the given directory contains an NMR experiment by verifying the existence of
a 'pdata/1/title' file, which is typical for Bruker NMR data structure.

Return `true` if the directory has the expected title file structure, `false` otherwise.
"""
function isexpt(directory)
    # check for existence of pdata/1/title file
    title_path = joinpath(directory, "pdata", "1", "title")
    return isfile(title_path)
end

function short_expt_path(directory)
    # shorten an experiment path to folder/number
    # e.g. /Users/chris/NMR/crick-701/sophia_trypsin_lig1_251117/10 -> sophia_trypsin_lig1_251117/10
    parent_folder = basename(dirname(directory))
    folder_num = basename(directory)
    return joinpath(parent_folder, folder_num)
end