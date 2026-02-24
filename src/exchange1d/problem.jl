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
    fit(prob::ExchangeProblem, params0::ComponentArray) -> NamedTuple

Fit all experiments jointly using least-squares optimisation.

Returns a NamedTuple with:
- `params`: fitted parameters with uncertainties (ComponentArray{Measurement})
- `params_value`: fitted parameters without uncertainties (ComponentArray{Float64}), for plotting
- `chi2`: chi-squared statistic
- `reduced_chi2`: chi-squared / degrees of freedom
- `covariance`: parameter covariance matrix
- `nobs`: number of observations
- `nparams`: number of fitted parameters
- `dof`: degrees of freedom
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
    covar = vcov(result)
    pfit_uncertain = ComponentArray(Measurements.correlated_values(result.param, covar), ax)

    # chi2 from weighted residuals
    predicted = model_func(dummy_x, result.param)
    chi2 = sum(((observed .- predicted) ./ errors) .^ 2)
    n_obs = length(observed)
    n_params = length(p0)
    dof = n_obs - n_params

    return (
        params = pfit_uncertain,
        params_value = pfit,
        chi2 = chi2,
        reduced_chi2 = chi2 / dof,
        covariance = covar,
        nobs = n_obs,
        nparams = n_params,
        dof = dof,
    )
end
