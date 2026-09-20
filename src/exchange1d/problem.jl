"""
    integrate!(prob::ExchangeProblem, peakppm, noiseppm, ppmwidth)

Integrate all experiments in the problem at the given peak and noise
positions, and record them on `prob.integration` for later reference (e.g.
when saving results).
"""
function integrate!(prob::ExchangeProblem, peakppm, noiseppm, ppmwidth)
    for expt in prob.experiments
        integrate!(expt, peakppm, noiseppm, ppmwidth)
    end
    prob.integration = (; peakppm, noiseppm, ppmwidth)
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
    isatbound(value, bound; tol=1e-6) -> Bool

Whether a fitted parameter has converged onto a finite bound rather than
settling to an interior optimum near it. `tol` is relative to the larger of
`value` and `bound` in magnitude (floored at 1).
"""
function isatbound(value::Float64, bound::Float64; tol::Float64=1e-6)
    isfinite(bound) || return false
    scale = max(abs(value), abs(bound), 1.0)
    return abs(value - bound) <= tol * scale
end

"""Correlation matrix corresponding to covariance matrix `cov`:
`cor[i,j] = cov[i,j] / sqrt(cov[i,i] * cov[j,j])`."""
function correlationmatrix(cov::Matrix{Float64})
    d = sqrt.(diag(cov))
    return cov ./ (d * d')
end

"""Off-diagonal `(i, j, r)` triples of `cor` with `|r| ≥ threshold` and
`i < j`, indexed into `cor`'s own row/column order (the fitted, non-fixed
parameters — see `FitResult.freeidx`), not the flat parameter indices used
elsewhere in `FitResult`."""
function strongcorrelations(cor::Matrix{Float64}; threshold::Float64=0.95)
    n = size(cor, 1)
    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:n, j in (i + 1):n
        r = cor[i, j]
        isfinite(r) && abs(r) >= threshold && push!(pairs, (i, j, r))
    end
    return pairs
end

"""
    warncorrelations(cor, freeitems, state_labels, fields)

Emit an `@warn` for each pair of fitted parameters correlated at or above
`strongcorrelations`'s threshold: this strong a correlation usually means the
two parameters trade off against each other rather than being separately
identifiable from this data, so their individually-reported uncertainties
understate that. `freeitems` are the `_ParamItem`s of the fitted (non-fixed)
parameters in `cor`'s row/column order. `_pretty_label` is defined in
interface.jl, included after this file — see `_ratelowerbound` above for why
that's fine.
"""
function warncorrelations(cor::Matrix{Float64}, freeitems, state_labels, fields)
    for (i, j, r) in strongcorrelations(cor)
        labeli = _pretty_label(freeitems[i], state_labels, fields)
        labelj = _pretty_label(freeitems[j], state_labels, fields)
        @warn "Strongly correlated fitted parameters ($labeli, $labelj): r = $(round(r; digits=3))"
    end
    return nothing
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

    items = _flatten_params_items(params0)
    lower = fill(-Inf, n)
    for item in items
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

    cor = correlationmatrix(covar)
    atbound = Set{Int}(freeidx[k] for k in eachindex(freeidx)
                       if isatbound(result.param[k], lower[freeidx[k]]))

    state_labels = states(prob.model)
    fields = _unique_fields(prob.experiments)
    warncorrelations(cor, items[freeidx], state_labels, fields)

    return FitResult(pfit_uncertain, pfit, ComponentArray(copy(p0), ax),
                     chi2, chi2 / dof, covar, cor,
                     n_obs, n_params, dof, copy(fixed), atbound, copy(freeidx), prob)
end
