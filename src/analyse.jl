"""
    AnalysisDispatch

A registry-based dispatch system for NMR analysis routines. This module provides:
- Classification of input files by types and features
- Rule-based matching of files to analysis routines
- Interactive multi-select menu for choosing analyses
- Extensible registration of new analysis rules
"""
module AnalysisDispatch

using REPL.TerminalMenus

export analyse, register_analysis!, MultiFileRule

"""
    ExperimentFile

Represents a classified NMR experiment file with its types and features.

# Fields
- `filename::String`: Path to the experiment file
- `types::Set{String}`: Set of experiment types (e.g., "1d", "calibration", "r1rho")
- `features::Set{String}`: Set of experiment features (e.g., "nutation", "on_resonance")
"""
struct ExperimentFile
    filename::String
    types::Set{String}
    features::Set{String}
end

function ExperimentFile(filename::String, types, features)
    return ExperimentFile(filename, Set{String}(types), Set{String}(features))
end

"""
    AnalysisRule

Abstract type for analysis rules that match experiments to analysis routines.
"""
abstract type AnalysisRule end

"""
    SingleFileRule <: AnalysisRule

A rule that matches individual files based on required types and features.

# Fields
- `required_types::Set{String}`: Types that must be present in the experiment
- `required_features::Set{String}`: Features that must be present in the experiment
- `handler::Function`: Function to call with matched file (receives ExperimentFile)
- `name::String`: Human-readable name for this analysis
"""
struct SingleFileRule <: AnalysisRule
    required_types::Set{String}
    required_features::Set{String}
    handler::Function
    name::String
end

function SingleFileRule(types, features, handler::Function, name::String)
    return SingleFileRule(Set{String}(types), Set{String}(features), handler, name)
end

"""
    MultiFileRule <: AnalysisRule

A rule that matches multiple files using a custom matcher function.

# Fields
- `matcher::Function`: Function that takes Vector{ExperimentFile} and returns
  matched files or nothing if no match
- `handler::Function`: Function to call with matched files (receives Vector{ExperimentFile})
- `name::String`: Human-readable name for this analysis
"""
struct MultiFileRule <: AnalysisRule
    matcher::Function
    handler::Function
    name::String
end

"""
    AnalysisOption

Represents a matched analysis that can be presented to the user.

# Fields
- `rule::AnalysisRule`: The rule that matched
- `matched_files::Vector{ExperimentFile}`: Files that matched this rule
- `description::String`: Human-readable description for menu display
"""
struct AnalysisOption
    rule::AnalysisRule
    matched_files::Vector{ExperimentFile}
    description::String
end

# Module-level registry
const ANALYSIS_REGISTRY = AnalysisRule[]

"""
    register_analysis!(rule::AnalysisRule)

Add an analysis rule to the registry.
"""
function register_analysis!(rule::AnalysisRule)
    push!(ANALYSIS_REGISTRY, rule)
    return rule
end

"""
    register_analysis!(required_types, required_features, handler, name)

Convenience function to register a SingleFileRule.

# Arguments
- `required_types`: Collection of required experiment types
- `required_features`: Collection of required experiment features
- `handler::Function`: Analysis function to call
- `name::String`: Human-readable name for the analysis

# Example
```julia
# Match files with types containing "1d" AND "calibration", no feature requirements
register_analysis!(["1d", "calibration"], [], analyse_1d_calibration, "1D calibration")

# Match files with type "r1rho" and feature "on_resonance"
register_analysis!(["r1rho"], ["on_resonance"], analyse_r1rho, "On-resonance R1ρ")
```
"""
function register_analysis!(required_types, required_features, handler::Function,
                            name::String)
    return register_analysis!(SingleFileRule(required_types, required_features, handler,
                                             name))
end

"""
    clear_registry!()

Clear all registered analysis rules. Primarily for testing.
"""
function clear_registry!()
    return empty!(ANALYSIS_REGISTRY)
end

"""
    match_rule(rule::SingleFileRule, experiments::Vector{ExperimentFile})

Match a single-file rule against experiments. Returns all experiments where
the required types and features are subsets of the experiment's types and features,
or `nothing` if no experiments match.
"""
function match_rule(rule::SingleFileRule, experiments::Vector{ExperimentFile})
    matched = filter(experiments) do exp
        return rule.required_types ⊆ exp.types && rule.required_features ⊆ exp.features
    end
    return isempty(matched) ? nothing : matched
end

"""
    match_rule(rule::MultiFileRule, experiments::Vector{ExperimentFile})

Match a multi-file rule against experiments by delegating to the rule's matcher function.
"""
function match_rule(rule::MultiFileRule, experiments::Vector{ExperimentFile})
    return rule.matcher(experiments)
end

"""
    describe(rule::SingleFileRule, file::ExperimentFile)

Generate a menu description for a single-file rule match.
"""
function describe(rule::SingleFileRule, file::ExperimentFile)
    return "$(rule.name): $(basename(file.filename))"
end

"""
    describe(rule::MultiFileRule, matched_files::Vector{ExperimentFile})

Generate a menu description for a multi-file rule match.
"""
function describe(rule::MultiFileRule, matched_files::Vector{ExperimentFile})
    return "$(rule.name): $(length(matched_files)) files"
end

"""
    specificity(option::AnalysisOption)

Calculate specificity score for sorting options.
Multi-file rules come first, then sorted by number of matched files descending.
"""
function specificity(option::AnalysisOption)
    is_multi = option.rule isa MultiFileRule
    n_files = length(option.matched_files)
    # Return tuple: multi-file rules first (true > false when negated), then by file count
    return (!is_multi, -n_files)
end

"""
    find_available_analyses(experiments::Vector{ExperimentFile})

Find all available analyses for the given experiments by matching against
the registry. Returns a sorted vector of AnalysisOption.

Single-file rules are expanded to one option per matching file.
Results are sorted by specificity: multi-file rules first, then by number
of matched files descending.
"""
function find_available_analyses(experiments::Vector{ExperimentFile})
    options = AnalysisOption[]

    for rule in ANALYSIS_REGISTRY
        matched = match_rule(rule, experiments)
        isnothing(matched) && continue

        if rule isa SingleFileRule
            # Expand to one option per matching file
            for file in matched
                desc = describe(rule, file)
                push!(options, AnalysisOption(rule, [file], desc))
            end
        else
            # Multi-file rule: single option with all matched files
            desc = describe(rule, matched)
            push!(options, AnalysisOption(rule, matched, desc))
        end
    end

    # Sort by specificity
    sort!(options; by=specificity)
    return options
end

"""
    run_analysis(option::AnalysisOption)

Execute the analysis for a given option.
"""
function run_analysis(option::AnalysisOption)
    if option.rule isa SingleFileRule
        # Single file rules receive the single matched file
        return option.rule.handler(option.matched_files[1])
    else
        # Multi file rules receive all matched files
        return option.rule.handler(option.matched_files)
    end
end

"""
    classify_file(filename::AbstractString, types_and_features_func::Function)

Classify a single file into an ExperimentFile using the provided classification function.
"""
function classify_file(filename::AbstractString, types_and_features_func::Function)
    types, features = types_and_features_func(filename)
    return ExperimentFile(filename, types, features)
end

"""
    analyse(filenames::AbstractVector{<:AbstractString}; types_and_features_func=Main.NMRAnalysis.types_and_features)

Analyse NMR experiment files using the registry-based dispatch system.

# Process
1. Classify files using `types_and_features_func`
2. Find available analyses from the registry
3. If no options, warn and return `nothing`
4. If one option, run it directly
5. Otherwise, present a multi-select menu and run selected analyses

# Arguments
- `filenames`: Vector of experiment file paths
- `types_and_features_func`: Function to classify files (default: NMRAnalysis.types_and_features)

# Returns
- Single result if one analysis was run
- Vector of results if multiple analyses were selected
- `nothing` if no analyses are available
"""
function analyse(filenames::AbstractVector{<:AbstractString};
                 types_and_features_func::Function=_default_types_and_features)
    # Classify all files
    experiments = [classify_file(f, types_and_features_func) for f in filenames]

    # Check for files that couldn't be classified
    unclassified = filter(e -> isempty(e.types), experiments)
    if !isempty(unclassified)
        for exp in unclassified
            @warn "Could not determine experiment type for file: $(exp.filename)"
        end
    end

    # Filter to only classified experiments
    classified = filter(e -> !isempty(e.types), experiments)
    if isempty(classified)
        @warn "No experiment types could be determined. Unable to proceed with analysis."
        return nothing
    end

    # Find available analyses
    options = find_available_analyses(classified)

    if isempty(options)
        types_found = unique(reduce(union, [e.types for e in classified]))
        @info "No automatic analysis routine available for experiment types: $(join(types_found, ", "))"
        return nothing
    end

    if length(options) == 1
        # Single option: run directly
        return run_analysis(options[1])
    end

    # Multiple options: present menu
    menu_items = [opt.description for opt in options]
    menu = MultiSelectMenu(menu_items)
    choices = request("Select analyses to run:", menu)

    if isempty(choices)
        @info "No analyses selected."
        return nothing
    end

    # Run selected analyses
    selected_indices = collect(choices)
    results = [run_analysis(options[i]) for i in selected_indices]

    return length(results) == 1 ? results[1] : results
end

"""
    analyse(filename::AbstractString; kwargs...)

Convenience wrapper for analysing a single file.
"""
function analyse(filename::AbstractString; kwargs...)
    return analyse([filename]; kwargs...)
end

# Get types_and_features from the parent NMRAnalysis module
function _default_types_and_features(filename)
    return parentmodule(@__MODULE__).types_and_features(filename)
end

end # module AnalysisDispatch

using .AnalysisDispatch
