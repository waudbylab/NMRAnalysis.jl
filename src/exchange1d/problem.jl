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
    residuals(expt::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})

Weighted residuals comparing the raw spin-lock decay intensities directly to the
Bloch-McConnell model, rather than to the rate independently fitted per condition into
`observed_intensities` (which exists only for display — see `plotresult!` and `ratewres`).

Each condition (one νSL on resonance, one offset off resonance) has its own equilibrium
intensity I₀. This is a linear parameter, separable from the nonlinear relaxation rate
that `simulate!` has just written into `predicted_intensities`, so rather than adding one
I₀ per condition to the fit, each is eliminated analytically by variable projection (Golub
& Pereyra, 1973, *SIAM J. Numer. Anal.* 10, 413–432): for a fixed rate R the weighted
least-squares I₀ minimising `Σ (I_obs - I₀ exp(-R t))²` has the closed form used below,
leaving only the exchange and relaxation parameters for the nonlinear optimiser.
"""
function residuals(expt::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})
    R = expt.predicted_intensities
    resid = similar(expt.rawintensities)
    for k in eachindex(R)
        y = @view expt.rawintensities[:, k]
        x = @. exp(-expt.TSL * R[k])
        I0 = dot(y, x) / sum(abs2, x)
        @. resid[:, k] = (y - I0 * x) / expt.rawnoise
    end
    return vec(resid)
end

"""
    nprofiledparams(expt::AbstractExperiment) -> Int

Number of parameters `residuals(expt)` fits implicitly by analytic profiling (variable
projection) rather than exposing them in the `ComponentArray` `fit` optimises. Zero by
default. R1ρ experiments override this to one I₀ per condition (see `residuals` above):
each is a genuine fitted parameter, so `fit` must count it towards `dof` and the parameter
covariance even though the optimiser never sees it directly.
"""
nprofiledparams(::AbstractExperiment) = 0
function nprofiledparams(expt::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})
    return size(expt.rawintensities, 2)
end

"""
    ratewres(expt::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})

Weighted residual of the rate independently fitted per condition (`observed_intensities`)
against the Bloch-McConnell prediction (`predicted_intensities`) — what the R1ρ result
plots actually draw. Distinct from `residuals(expt)`, which compares raw intensities
directly and is what the joint fit minimises.
"""
function ratewres(expt::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})
    yobs = expt.observed_intensities
    return (Measurements.value.(yobs) .- expt.predicted_intensities) ./
           Measurements.uncertainty.(yobs)
end

"""
    residuals(prob::ExchangeProblem, params::ComponentArray)

Simulate all experiments then return concatenated weighted residuals. This is the
objective the joint fit minimises, so a per-experiment `residuals` override (e.g. R1ρ's,
above) changes what the fit itself compares, not just what is reported afterwards.
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

The objective minimised is `residuals(prob, params)` — the concatenation of each
experiment's own weighted residuals — evaluated against a zero target, so a
per-experiment `residuals` override (as R1ρ uses to eliminate I₀ by variable projection;
see `residuals(::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})`) changes what the
fit itself minimises, not just what is reported afterwards.
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

    function paramsfor(pfree)
        pfull = copy(p0)
        pfull[freeidx] .= pfree
        return ComponentArray(pfull, ax)
    end
    resid_func(x, pfree) = residuals(prob, paramsfor(pfree))

    n_obs = length(resid_func(nothing, p0[freeidx]))
    dummy_x = 1:n_obs
    target = zeros(n_obs)

    result = curve_fit(resid_func, dummy_x, target, p0[freeidx]; lower=lower[freeidx])

    # R1ρ profiles one I₀ per condition out of residuals(expt) by variable projection (see
    # residuals(::Union{R1rhoOnResExperiment,R1rhoOffResExperiment})) rather than exposing
    # it to curve_fit, but each is still a genuine fitted parameter and must be counted
    # here — otherwise dof, reduced_chi2 and (below) the parameter covariance are all
    # biased as if those amplitudes came for free
    n_profiled = sum(nprofiledparams(expt) for expt in prob.experiments)
    n_params = length(freeidx) + n_profiled
    dof = n_obs - n_params

    # reconstruct as ComponentArrays, re-inserting fixed values
    pfull_fit = copy(p0)
    pfull_fit[freeidx] .= result.param
    pfit = ComponentArray(pfull_fit, ax)

    covar = try
        # curve_fit's own vcov scales by chi2 / (n_obs - length(freeidx)), oblivious to the
        # profiled parameters above; its inv(J'J) shape is unaffected (that's the useful
        # property of variable projection) but the noise-variance scalar needs rescaling
        # onto the true dof for the reported uncertainties to reflect every fitted parameter
        vcov(result) .* (n_obs - length(freeidx)) / dof
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

    chi2 = sum(result.resid .^ 2)

    d = sqrt.(diag(covar))
    cor = covar ./ (d * d')
    atbound = Set{Int}(freeidx[k]
                       for k in eachindex(freeidx)
                       if isatbound(result.param[k], lower[freeidx[k]]))

    state_labels = states(prob.model)
    fields = _unique_fields(prob.experiments)
    warncorrelations(cor, items[freeidx], state_labels, fields)

    return FitResult(pfit_uncertain, pfit, ComponentArray(copy(p0), ax),
                     chi2, chi2 / dof, covar, cor,
                     n_obs, n_params, dof, copy(fixed), atbound, copy(freeidx), prob)
end
