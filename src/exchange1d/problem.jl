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
    _ratelowerbound(item) -> Float64

Lower bound for a raw (linear, un-fixed) parameter passed to the optimiser.
Relaxation rates (`spin.R1_*`, `spin.R2_*`) are physically non-negative but,
unlike Kd/koff (see `_islogparam`) or populations (see `_isdgparam`), are
stored and fitted directly in linear space — so unlike those, nothing stops
the optimiser from stepping into negative territory mid-search, which makes
the Bloch-McConnell matrix exponential used to simulate CEST/R1ρ experiments
numerically unstable (a relaxation rate < 0 is a growing, not decaying,
mode). Bounding these at zero keeps the search in the physical region
without changing how they're entered, displayed, or stored. Takes an
`_ParamItem` (from `_flatten_params_items`, defined in interface.jl); left
untyped here since interface.jl is included after this file.
"""
function _ratelowerbound(item)
    item.section == "spin" || return -Inf
    return match(r"^R[12]_", String(_paramkey(item))) !== nothing ? 0.0 : -Inf
end

"""
    fit(prob::ExchangeProblem, params0::ComponentArray; fixed=Set{Int}()) -> FitResult

Fit all experiments jointly using least-squares optimisation.

`fixed` is a set of flat indices (as produced by `_flatten_params_items`) into
`params0` that are held constant at their `params0` value rather than being
optimised. This lets, e.g., a chemical shift be pinned by the user while the
rest of the model is fitted.

Returns a `FitResult` containing fitted parameters (with uncertainties),
fit statistics, and a reference to the problem for display and plotting.
"""
function fit(prob::ExchangeProblem, params0::ComponentArray; fixed::Set{Int}=Set{Int}())
    p0 = collect(params0)
    ax = getaxes(params0)
    n = length(p0)
    freeidx = [i for i in 1:n if i ∉ fixed]

    lower = fill(-Inf, n)
    for item in _flatten_params_items(params0)
        lower[item.flat_index] = _ratelowerbound(item)
    end

    # observed values and weights from all experiments
    observed = vcat([Measurements.value.(expt.observed_intensities)
                     for expt in prob.experiments]...)
    errors = vcat([Measurements.uncertainty.(expt.observed_intensities)
                   for expt in prob.experiments]...)
    wt = errors .^ -2

    dummy_x = 1:length(observed)

    function model_func(x, pfree)
        pfull = copy(p0)
        pfull[freeidx] .= pfree
        params = ComponentArray(pfull, ax)
        simulate!(prob, params)
        return vcat([copy(expt.predicted_intensities) for expt in prob.experiments]...)
    end

    result = curve_fit(model_func, dummy_x, observed, wt, p0[freeidx]; lower=lower[freeidx])

    # reconstruct as ComponentArrays, re-inserting fixed values
    pfull_fit = copy(p0)
    pfull_fit[freeidx] .= result.param
    pfit = ComponentArray(pfull_fit, ax)

    covar = try
        vcov(result)
    catch e
        @error "Failed to compute covariance matrix"
        zeros(length(freeidx), length(freeidx))
    end

    full_uncertain = Vector{Measurement{Float64}}(undef, n)
    full_uncertain[freeidx] .= Measurements.correlated_values(result.param, covar)
    for i in fixed
        full_uncertain[i] = p0[i] ± 0.0
    end
    pfit_uncertain = ComponentArray(full_uncertain, ax)

    # chi2 from weighted residuals
    predicted = model_func(dummy_x, result.param)
    chi2 = sum(((observed .- predicted) ./ errors) .^ 2)
    n_obs = length(observed)
    n_params = length(freeidx)
    dof = n_obs - n_params

    return FitResult(pfit_uncertain, pfit, ComponentArray(copy(p0), ax),
                     chi2, chi2 / dof, covar,
                     n_obs, n_params, dof, copy(fixed), prob)
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
