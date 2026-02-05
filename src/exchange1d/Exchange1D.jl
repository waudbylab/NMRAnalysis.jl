"""
    Exchange1D

Module for 1D chemical exchange analysis using Bloch-McConnell equations.

Supports R1ρ relaxation dispersion and CEST experiments with:
- NoExchange (null model) and TwoState exchange models
- Full Bloch-McConnell simulation
- Model comparison via F-test and AIC
"""
module Exchange1D

using Distributions: cdf, FDist
using LinearAlgebra
using LsqFit
using Measurements: Measurements
using NMRTools
using Plots
using REPL.TerminalMenus
using Statistics

export exchange1d

# Import from parent module
using ..NMRAnalysis: analyse, register_analysis!, MultiFileRule
using ..NMRAnalysis: get1dregionandnoise, relaxation1d

# Include submodules in dependency order
include("parameters.jl")
include("models.jl")
include("bloch_mcconnell.jl")
include("experiments.jl")
include("predict.jl")
include("fitting.jl")
include("loading.jl")
include("plotting.jl")

# Export types and key functions
export Parameter, fixed, value, LOG_TRANSFORM
export AbstractModel, NoExchange, TwoState, nstates, modelname
export AbstractExperiment, R1rhoExperiment, CESTExperiment
export predict, fit_exchange, compare_models

"""
    exchange1d(filenames::Vector{String})

Analyze 1D chemical exchange experiments (R1ρ and/or CEST).

# Workflow
1. Load and classify experiments from provided files
2. Interactively select integration and noise regions
3. Integrate peaks to get intensities
4. Use R1 calibration if available
5. Select exchange models to fit (multi-select)
6. Fit each model with user-provided initial parameters
7. Compare models and display results

# Arguments
- `filenames`: Vector of paths to NMR experiment directories

# Returns
Named tuple with:
- `results`: Dict mapping models to FitResult
- `plt`: Summary plot
"""
function exchange1d(filenames::Vector{String})
    @info "Starting exchange analysis..."

    # 1. Classify experiments
    r1rho_files = filter_r1rho_files(filenames)
    cest_files = filter_cest_files(filenames)
    r1_files = filter_r1_calibration_files(filenames)

    @info "Found $(length(r1rho_files)) R1ρ, $(length(cest_files)) CEST, $(length(r1_files)) R1 calibration experiments"

    if isempty(r1rho_files) && isempty(cest_files)
        error("No R1ρ or CEST experiments found in provided files")
    end

    # 2. Load data
    r1rho_data = isempty(r1rho_files) ? nothing : load_r1rho_data(r1rho_files)
    cest_data = isempty(cest_files) ? nothing : load_cest_data(cest_files)

    # 3. Interactive region selection
    slice = if !isnothing(r1rho_data)
        r1rho_data.spectra[1]
    else
        cest_data.spectra[1]
    end
    roi, noiseroi = get1dregionandnoise(slice)

    # 4. Integrate all experiments
    experiments = AbstractExperiment[]
    if !isnothing(r1rho_data)
        intensities, uncertainties = integrate_spectra(r1rho_data, roi, noiseroi)
        push!(experiments, make_r1rho_experiment(r1rho_data, intensities, uncertainties))
        @info "Integrated $(length(intensities)) R1ρ data points"
    end
    if !isnothing(cest_data)
        intensities, uncertainties = integrate_spectra(cest_data, roi, noiseroi)
        push!(experiments, make_cest_experiment(cest_data, intensities, uncertainties))
        @info "Integrated $(length(intensities)) CEST data points"
    end

    # 5. R1 calibration (if available)
    R1 = if !isempty(r1_files)
        r1_result = relaxation1d(r1_files[1], roi, noiseroi)
        R1_val = Measurements.value(r1_result.rate)
        @info "Using R1 = $(round(R1_val, digits=3)) s⁻¹ from calibration"
        R1_val
    else
        nothing
    end

    # 6. Model selection via terminal menu (multi-select)
    models = select_models_menu()

    # 7. Fit and compare models
    results = fit_and_compare_models(experiments, models, R1)

    return results
end

"""
    select_models_menu()

Display multi-select menu for model selection.
"""
function select_models_menu()
    options = ["No exchange (null model)", "Two-state exchange"]
    menu = MultiSelectMenu(options)
    choices = request("Select models to fit (space to select, enter to confirm):", menu)

    models = AbstractModel[]
    1 in choices && push!(models, NoExchange())
    2 in choices && push!(models, TwoState())

    if isempty(models)
        error("No models selected")
    end

    return models
end

"""
    fit_and_compare_models(experiments, models::Vector{<:AbstractModel}, R1)

Fit each selected model, prompting user for initial parameters, then compare results.
"""
function fit_and_compare_models(experiments::Vector{<:AbstractExperiment},
                                 models::Vector{<:AbstractModel}, R1)
    results = Dict{AbstractModel, FitResult}()

    for model in models
        println("\n" * "="^50)
        @info "Setting up $(modelname(model)) model..."

        # Get model-specific parameters from user
        modelpars = prompt_model_parameters(model)

        # Build spin parameters (common across models, but nstates differs)
        spinpars = prompt_spin_parameters(model, experiments, R1)

        # Fit
        @info "Fitting $(modelname(model))..."
        results[model] = fit_exchange(experiments, model, modelpars, spinpars)

        # Display results for this model
        display_fit_result(results[model])
    end

    # Compare models if more than one fitted
    if length(models) > 1
        comparison = compare_models(results)
        display_comparison(comparison)
    end

    # Plot results
    plt = plot_results(experiments, results)
    display(plt)

    return (results=results, plt=plt)
end

"""
    prompt_model_parameters(model::NoExchange)

No parameters needed for NoExchange model.
"""
function prompt_model_parameters(::NoExchange)
    NamedTuple()
end

"""
    prompt_model_parameters(model::TwoState)

Prompt user for two-state exchange parameters (kex, pB).
"""
function prompt_model_parameters(::TwoState)
    println("Enter initial parameter estimates for two-state exchange:")

    print("  Exchange rate kex (s⁻¹) [default 500]: ")
    kex_input = readline()
    kex = isempty(strip(kex_input)) ? 500.0 : parse(Float64, kex_input)

    print("  Minor state population pB [default 0.05]: ")
    pB_input = readline()
    pB = isempty(strip(pB_input)) ? 0.05 : parse(Float64, pB_input)

    (
        kex = Parameter(kex; transform=LOG_TRANSFORM, bounds=(1.0, 1e6)),
        pB = Parameter(pB; bounds=(0.001, 0.5)),
    )
end

"""
    prompt_spin_parameters(model::AbstractModel, experiments, R1)

Prompt user for spin parameters (R1, R2, δ for each state).
"""
function prompt_spin_parameters(model::AbstractModel,
                                 experiments::Vector{<:AbstractExperiment}, R1)
    n = nstates(model)
    println("Enter spin parameters:")

    # R1 - use calibration if available, otherwise prompt
    R1_val = if !isnothing(R1)
        println("  R1 = $(round(R1, digits=3)) s⁻¹ (from calibration, fixed)")
        R1
    else
        print("  R1 (s⁻¹) [default 1.5]: ")
        input = readline()
        isempty(strip(input)) ? 1.5 : parse(Float64, input)
    end

    # R2 and δ for each state
    R2_params = Parameter[]
    δ_params = Parameter[]

    for i in 1:n
        state_name = n == 1 ? "" : (i == 1 ? " (ground state)" : " (excited state)")

        print("  R2$state_name (s⁻¹) [default 15]: ")
        R2_input = readline()
        R2_val = isempty(strip(R2_input)) ? 15.0 : parse(Float64, R2_input)
        push!(R2_params, Parameter(R2_val; bounds=(0.0, 1000.0)))

        print("  Chemical shift δ$state_name (rad/s) [default 0]: ")
        δ_input = readline()
        δ_val = isempty(strip(δ_input)) ? 0.0 : parse(Float64, δ_input)
        push!(δ_params, Parameter(δ_val))
    end

    # Estimate amplitude from data
    amp_estimate = estimate_amplitude(experiments)

    (
        R1 = fixed(R1_val),
        R2 = Tuple(R2_params),
        δ = Tuple(δ_params),
        R1rho_amplitude = Parameter(amp_estimate; transform=LOG_TRANSFORM),
        CEST_amplitude = Parameter(amp_estimate; transform=LOG_TRANSFORM),
    )
end

"""
    estimate_amplitude(experiments)

Estimate initial amplitude from maximum observed intensity.
"""
function estimate_amplitude(experiments::Vector{<:AbstractExperiment})
    max_intensity = 0.0
    for exp in experiments
        if hasdata(exp)
            max_intensity = max(max_intensity, maximum(exp.intensity))
        end
    end
    max_intensity > 0 ? max_intensity : 1.0
end

# Registration with analysis system
function __init__()
    rule = MultiFileRule(
        expts -> begin
            oneD = filter(e -> "1d" in e.types, expts)
            r1rho = filter(e -> "r1rho" in e.types, oneD)
            cest = filter(e -> "cest" in e.types, oneD)
            r1cal = filter(e -> "relaxation" in e.types && "R1" in e.features, oneD)
            combined = vcat(r1rho, cest, r1cal)
            # Only match if we have R1ρ or CEST experiments
            (length(r1rho) > 0 || length(cest) > 0) ? combined : nothing
        end,
        expts -> exchange1d([e.filename for e in expts]),
        "Exchange analysis (R1ρ and CEST)"
    )
    return register_analysis!(rule)
end

end # module
