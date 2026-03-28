"""
    integrate!(prob::ExchangeProblem, peakppm, noiseppm, ppmwidth)

Integrate all experiments in the problem at the given peak and noise positions.
"""
function integrate!(prob::ExchangeProblem, peakppm, noiseppm, ppmwidth)
    for expt in prob.experiments
        integrate!(expt, peakppm, noiseppm, ppmwidth)
    end
    return nothing
end

"""
    simulate!(prob::ExchangeProblem, params::ComponentArray)

Simulate predicted values for all experiments in the problem.
"""
function simulate!(prob::ExchangeProblem, params::ComponentArray)
    for expt in prob.experiments
        simulate!(expt, prob.model, params)
    end
    return nothing
end

"""
    residuals(expt::AbstractExperiment)

Return weighted residuals `(observed - predicted) / uncertainty` for an experiment.
Default implementation using the `observed_intensities` and `predicted_intensities` fields.
"""
function residuals(expt::AbstractExperiment)
    obs_values = Measurements.value.(expt.observed_intensities)
    obs_errors = Measurements.uncertainty.(expt.observed_intensities)
    return (obs_values .- expt.predicted_intensities) ./ obs_errors
end

"""
    residuals(prob::ExchangeProblem, params::ComponentArray)

Simulate all experiments then return concatenated weighted residuals.
"""
function residuals(prob::ExchangeProblem, params::ComponentArray)
    simulate!(prob, params)
    return vcat([residuals(expt) for expt in prob.experiments]...)
end

"""
    fit(prob::ExchangeProblem, params0::ComponentArray) -> FitResult

Fit all experiments jointly using least-squares optimisation.

Returns a `FitResult` containing fitted parameters (with uncertainties),
fit statistics, and a reference to the problem for display and plotting.
"""
function fit(prob::ExchangeProblem, params0::ComponentArray)
    p0 = collect(params0)
    ax = getaxes(params0)

    # observed values and weights from all experiments
    observed = vcat([Measurements.value.(expt.observed_intensities)
                     for expt in prob.experiments]...)
    errors = vcat([Measurements.uncertainty.(expt.observed_intensities)
                   for expt in prob.experiments]...)
    wt = errors .^ -2

    dummy_x = 1:length(observed)

    function model_func(x, p)
        params = ComponentArray(p, ax)
        simulate!(prob, params)
        return vcat([copy(expt.predicted_intensities) for expt in prob.experiments]...)
    end

    result = curve_fit(model_func, dummy_x, observed, wt, p0)

    # reconstruct as ComponentArrays
    pfit = ComponentArray(result.param, ax)
    covar = try
        vcov(result)
    catch e
        @error "Failed to compute covariance matrix"
        zeros(length(p0), length(p0))
    end
    pfit_uncertain = ComponentArray(Measurements.correlated_values(result.param, covar), ax)

    # chi2 from weighted residuals
    predicted = model_func(dummy_x, result.param)
    chi2 = sum(((observed .- predicted) ./ errors) .^ 2)
    n_obs = length(observed)
    n_params = length(p0)
    dof = n_obs - n_params

    return FitResult(pfit_uncertain, pfit, ComponentArray(copy(p0), ax),
                     chi2, chi2 / dof, covar,
                     n_obs, n_params, dof, prob)
end

"""
    plot_result(prob::ExchangeProblem, fit_result)

Plot all experiments in the problem using the fitted parameters from `fit_result`
(as returned by `fit`). Returns a vector of plots, one per experiment.
"""
function plot_result(prob::ExchangeProblem, fit_result, kwargs...)
    params = fit_result.params_value
    simulate!(prob, params)  # update predicted_intensities for all experiments
    plots = [plot_result(expt, fit_result; kwargs...) for expt in prob.experiments]
    return plots
end
