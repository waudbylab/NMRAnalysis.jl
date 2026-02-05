# Fitting infrastructure for exchange analysis

using LsqFit
using Distributions: cdf, FDist

"""
    PackedParameters

Metadata for unpacking fitted parameters back to NamedTuples.
"""
struct PackedParameters
    indices::Vector{Tuple{Symbol, Any}}  # (source, key) for each packed parameter
    n_model::Int                          # Number of model parameters
    n_spin::Int                           # Number of spin parameters
end

"""
    pack(modelpars, spinpars)

Pack non-fixed parameters into a flat vector for optimization.

Parameters are transformed to fitting space using their transform functions.

# Returns
- `x0`: Vector of initial parameter values in fitting space
- `metadata`: PackedParameters for unpacking
"""
function pack(modelpars, spinpars)
    params = Float64[]
    indices = Tuple{Symbol, Any}[]

    # Pack model parameters
    for (name, p) in pairs(modelpars)
        if !p.fixed
            push!(params, to_fit(p))
            push!(indices, (:model, name))
        end
    end
    n_model = length(params)

    # Pack spin parameters (handle nested tuples)
    for (name, val) in pairs(spinpars)
        if val isa Parameter
            if !val.fixed
                push!(params, to_fit(val))
                push!(indices, (:spin, name))
            end
        elseif val isa Tuple  # e.g., δ = (Parameter, Parameter)
            for (i, p) in enumerate(val)
                if !p.fixed
                    push!(params, to_fit(p))
                    push!(indices, (:spin, (name, i)))
                end
            end
        end
    end
    n_spin = length(params) - n_model

    return params, PackedParameters(indices, n_model, n_spin)
end

"""
    unpack(x, metadata, modelpars, spinpars)

Unpack parameter vector back to NamedTuples with fitted values.

# Returns
- `modelpars_fitted`: NamedTuple with fitted model parameters
- `spinpars_fitted`: NamedTuple with fitted spin parameters
"""
function unpack(x, metadata::PackedParameters, modelpars, spinpars)
    # Create mutable copies as Dicts
    mp = Dict{Symbol, Parameter}(pairs(modelpars))
    sp = Dict{Symbol, Any}(pairs(spinpars))

    for (i, (src, key)) in enumerate(metadata.indices)
        if src == :model
            p = mp[key]
            fitted_val = from_fit(p, x[i])
            mp[key] = with_fitted(p, fitted_val)
        else  # :spin
            if key isa Symbol
                p = sp[key]
                fitted_val = from_fit(p, x[i])
                sp[key] = with_fitted(p, fitted_val)
            else  # Tuple index (name, idx)
                name, idx = key
                tup = collect(sp[name])
                p = tup[idx]
                fitted_val = from_fit(p, x[i])
                tup[idx] = with_fitted(p, fitted_val)
                sp[name] = Tuple(tup)
            end
        end
    end

    return NamedTuple(mp), NamedTuple(sp)
end

"""
    get_bounds(modelpars, spinpars, metadata)

Extract lower and upper bounds for packed parameters.
"""
function get_bounds(modelpars, spinpars, metadata::PackedParameters)
    lower = Float64[]
    upper = Float64[]

    for (src, key) in metadata.indices
        if src == :model
            p = modelpars[key]
        else
            if key isa Symbol
                p = spinpars[key]
            else
                name, idx = key
                p = spinpars[name][idx]
            end
        end

        # Transform bounds to fitting space
        lb, ub = p.bounds

        # Handle edge cases for transforms
        # For log transform, 0 and negative values map to -Inf, Inf maps to Inf
        lb_trans = try
            p.transform[1](lb)
        catch
            -Inf
        end
        ub_trans = try
            p.transform[1](ub)
        catch
            Inf
        end

        push!(lower, lb_trans)
        push!(upper, ub_trans)
    end

    return lower, upper
end

"""
    residuals(exp::AbstractExperiment, model::AbstractModel, modelpars, spinpars)

Compute weighted residuals for an experiment.

Returns empty vector if experiment has no observed data.
"""
function residuals(exp::AbstractExperiment, model::AbstractModel, modelpars, spinpars)
    !hasdata(exp) && return Float64[]

    pred = predict(exp, model, modelpars, spinpars)
    (exp.intensity .- pred) ./ exp.σ
end

"""
    FitResult

Results from fitting an exchange model.

# Fields
- `model`: The fitted model
- `modelpars`: Fitted model parameters
- `spinpars`: Fitted spin parameters
- `χ2`: Sum of squared residuals
- `ndata`: Number of data points
- `npars`: Number of fitted parameters
- `predictions`: Predicted intensities for each experiment
"""
struct FitResult
    model::AbstractModel
    modelpars::NamedTuple
    spinpars::NamedTuple
    χ2::Float64
    ndata::Int
    npars::Int
    predictions::Vector{Vector{Float64}}
end

"""
    fit_exchange(experiments, model::AbstractModel, modelpars, spinpars)

Fit an exchange model to experimental data.

# Arguments
- `experiments`: Vector of experiments (R1rhoExperiment or CESTExperiment)
- `model`: Exchange model to fit
- `modelpars`: Initial model parameters
- `spinpars`: Initial spin parameters

# Returns
- `FitResult` with fitted parameters and statistics
"""
function fit_exchange(experiments::Vector{<:AbstractExperiment}, model::AbstractModel,
                      modelpars, spinpars)
    # Pack parameters
    x0, metadata = pack(modelpars, spinpars)

    # Count data points
    ndata = sum(npoints(exp) for exp in experiments if hasdata(exp))

    # If no free parameters, just compute predictions
    if isempty(x0)
        predictions = [predict(exp, model, modelpars, spinpars) for exp in experiments]
        r = vcat([residuals(exp, model, modelpars, spinpars) for exp in experiments]...)
        χ2 = sum(abs2, r)
        return FitResult(model, modelpars, spinpars, χ2, ndata, 0, predictions)
    end

    # Objective function
    function objective(x)
        mp, sp = unpack(x, metadata, modelpars, spinpars)
        r = vcat([residuals(exp, model, mp, sp) for exp in experiments]...)
        sum(abs2, r)
    end

    # Get bounds
    lower, upper = get_bounds(modelpars, spinpars, metadata)

    # Build combined intensity and uncertainty vectors for LsqFit
    ydata = vcat([exp.intensity for exp in experiments if hasdata(exp)]...)
    weights = vcat([1.0 ./ exp.σ for exp in experiments if hasdata(exp)]...)

    # Model function for LsqFit (dummy x argument)
    function model_func(_, p)
        mp, sp = unpack(p, metadata, modelpars, spinpars)
        vcat([predict(exp, model, mp, sp) for exp in experiments if hasdata(exp)]...)
    end

    # Fit using LsqFit
    fit = curve_fit(model_func, 1:ndata, ydata, weights.^2, x0; lower=lower, upper=upper)

    if !fit.converged
        @warn "Fit did not converge"
    end

    # Unpack results
    mp_fit, sp_fit = unpack(coef(fit), metadata, modelpars, spinpars)

    # Compute predictions and χ²
    predictions = [predict(exp, model, mp_fit, sp_fit) for exp in experiments]
    r = vcat([residuals(exp, model, mp_fit, sp_fit) for exp in experiments]...)
    χ2 = sum(abs2, r)

    return FitResult(model, mp_fit, sp_fit, χ2, ndata, length(x0), predictions)
end

"""
    compare_models(results::Dict{<:AbstractModel, FitResult})

Compare fitted models using F-test and AIC.

# Returns
Named tuple with:
- `comparisons`: Vector of pairwise F-test results
- `aic`: Dict of AIC values for each model
"""
function compare_models(results::Dict{<:AbstractModel, FitResult})
    comparisons = NamedTuple[]

    # Pairwise F-test for nested models (NoExchange nested in TwoState)
    no_exchange_key = findfirst(k -> k isa NoExchange, keys(results))
    two_state_key = findfirst(k -> k isa TwoState, keys(results))

    if !isnothing(no_exchange_key) && !isnothing(two_state_key)
        null_result = results[no_exchange_key]
        full_result = results[two_state_key]

        rss_null = null_result.χ2
        rss_full = full_result.χ2
        ndata = full_result.ndata
        df_null = ndata - null_result.npars
        df_full = ndata - full_result.npars

        # F-statistic
        F = ((rss_null - rss_full) / (df_null - df_full)) / (rss_full / df_full)
        p_value = 1 - cdf(FDist(df_null - df_full, df_full), F)

        push!(comparisons, (
            null_model = NoExchange(),
            full_model = TwoState(),
            F = F,
            p_value = p_value,
        ))
    end

    # AIC for all models
    aic_values = Dict{AbstractModel, Float64}()
    for (model, result) in results
        ndata = result.ndata
        # AIC = n*ln(RSS/n) + 2k
        aic = ndata * log(result.χ2 / ndata) + 2 * result.npars
        aic_values[model] = aic
    end

    return (comparisons=comparisons, aic=aic_values)
end

"""
    display_fit_result(result::FitResult)

Display fitted parameters and fit statistics.
"""
function display_fit_result(result::FitResult)
    println("\n=== $(modelname(result.model)) Fit Results ===")

    # Model parameters
    if !isempty(result.modelpars)
        println("Model parameters:")
        for (name, p) in pairs(result.modelpars)
            println("  $name = $(round(value(p), sigdigits=4))")
        end
    end

    # Key spin parameters
    println("Spin parameters:")
    println("  R1 = $(round(value(result.spinpars.R1), sigdigits=4)) s⁻¹")
    for (i, p) in enumerate(result.spinpars.R2)
        state = length(result.spinpars.R2) == 1 ? "" : " (state $i)"
        println("  R2$state = $(round(value(p), sigdigits=4)) s⁻¹")
    end
    for (i, p) in enumerate(result.spinpars.δ)
        state = length(result.spinpars.δ) == 1 ? "" : " (state $i)"
        println("  δ$state = $(round(value(p), sigdigits=4)) rad/s")
    end

    # Fit statistics
    println("Fit statistics:")
    println("  χ² = $(round(result.χ2, digits=2))")
    println("  Reduced χ² = $(round(result.χ2 / (result.ndata - result.npars), digits=3))")
    println("  Data points = $(result.ndata)")
    println("  Free parameters = $(result.npars)")
end

"""
    display_comparison(comparison)

Display model comparison results.
"""
function display_comparison(comparison)
    println("\n=== Model Comparison ===")

    # F-test results
    for c in comparison.comparisons
        println("\nF-test: $(modelname(c.null_model)) vs $(modelname(c.full_model))")
        println("  F = $(round(c.F, digits=2))")
        println("  p-value = $(round(c.p_value, sigdigits=3))")
        if c.p_value < 0.05
            println("  → Exchange is statistically significant (p < 0.05)")
        else
            println("  → No significant evidence for exchange (p ≥ 0.05)")
        end
    end

    # AIC comparison
    println("\nAIC values:")
    for (model, aic) in comparison.aic
        println("  $(modelname(model)): $(round(aic, digits=1))")
    end
    if length(comparison.aic) > 1
        best_model, best_aic = findmin(comparison.aic)
        println("  → Best model by AIC: $(modelname(best_model))")
    end
end
