# Analysis Rules

The `analyse()` function uses a registry-based dispatch system to match experiment files to appropriate analysis routines. This page describes how to extend the system with custom analysis rules.

## How It Works

1. **Classification**: Each input file is classified by its `types` (e.g., "1d", "calibration", "r1rho") and `features` (e.g., "nutation", "on_resonance")
2. **Matching**: Registered rules are matched against the classified files
3. **Selection**: If multiple analyses match, an interactive menu is presented
4. **Execution**: Selected analyses are run and results returned

## Registering Single-File Rules

For analyses that operate on individual files, use the convenience function:

```julia
register_analysis!(required_types, required_features, handler, name)
```

**Arguments:**
- `required_types`: Types that must be present in the experiment (matched as subset)
- `required_features`: Features that must be present (matched as subset)
- `handler`: Function receiving an `ExperimentFile` object
- `name`: Human-readable name shown in menus

**Example:**
```julia
# Match files where types include both "1d" and "calibration"
register_analysis!(
    ["1d", "calibration"],
    [],
    exp -> analyse_1d_calibration(exp.filename),
    "1D calibration"
)

# Match files with type "r1rho" and feature "on_resonance"
register_analysis!(
    ["r1rho"],
    ["on_resonance"],
    exp -> analyse_onres_r1rho(exp.filename),
    "On-resonance R1ρ"
)
```

The handler receives an `ExperimentFile` with fields:
- `filename::String`: Path to the experiment
- `types::Set{String}`: Experiment types
- `features::Set{String}`: Experiment features

## Registering Multi-File Rules

For analyses requiring multiple files (e.g., dispersion experiments), use `MultiFileRule`:

```julia
register_analysis!(MultiFileRule(matcher, handler, name))
```

**Arguments:**
- `matcher`: Function taking `Vector{ExperimentFile}`, returning matched files or `nothing`
- `handler`: Function receiving `Vector{ExperimentFile}`
- `name`: Human-readable name

**Example:**
```julia
register_analysis!(MultiFileRule(
    # Matcher: require at least 3 on-resonance R1ρ experiments
    experiments -> begin
        matched = filter(experiments) do e
            "r1rho" in e.types && "on_resonance" in e.features
        end
        length(matched) >= 3 ? matched : nothing
    end,
    # Handler: analyse all matched files together
    files -> analyse_r1rho_dispersion([f.filename for f in files]),
    "On-resonance R1ρ dispersion"
))
```

## Menu Ordering

Analysis options are sorted by specificity:
1. Multi-file rules appear first
2. Then sorted by number of matched files (descending)

This ensures more specific analyses (requiring multiple files) are prioritized over generic single-file analyses.

## Example: Complete Registration

```julia
using NMRAnalysis

# Single-file: basic 1D processing
register_analysis!(["1d"], [], process_1d, "Process 1D")

# Single-file: specific calibration analysis
register_analysis!(["1d", "calibration"], ["nutation"], analyse_nutation, "Nutation calibration")

# Multi-file: relaxation series
register_analysis!(MultiFileRule(
    experiments -> begin
        matched = filter(e -> "relaxation" in e.types, experiments)
        length(matched) >= 2 ? matched : nothing
    end,
    analyse_relaxation_series,
    "Relaxation series"
))
```

When `analyse()` is called with files matching multiple rules, users see a menu like:

```
Select analyses to run:
 [ ] Relaxation series: 5 files
 [ ] Nutation calibration: expt_001
 [ ] Process 1D: expt_001
 [ ] Process 1D: expt_002
```
