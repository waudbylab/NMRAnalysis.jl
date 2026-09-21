abstract type AbstractModel end

"""
    modelorder(model::AbstractModel) -> Int

Sort key controlling where `model` appears in the model selection menu
(lowest first). Defaults to placing a model without an explicit ordering
after every model that has one, so a custom model added by a package
extension (see `docs/src/advanced/extending_exchange1d.md`) still appears
automatically without needing to touch this function — just define
`modelorder` too if a particular position matters.
"""
modelorder(::AbstractModel) = typemax(Int)

abstract type AbstractExperiment end

"""
    ExchangeProblem

A set of experiments to be fitted jointly to a shared exchange `model`.

`integration` records the peak-picking parameters (`peakppm`, `noiseppm`,
`ppmwidth`) used to integrate all experiments in the problem, as a
`NamedTuple`, or `nothing` before `integrate!` has been called. Kept here
(rather than discarded once used) so it can be saved alongside the fit
results for traceability back to the source spectra.

`savecalibration` is how a B₁ calibration fitted for this problem is kept:
a function writing that analysis into a folder, or `nothing` where the
calibration was given ready-made or not at all. Saving it is deferred to
`_save_results`, so the calibration lands inside the output folder of the
fit that used it, and a fit that is never saved leaves nothing behind.
"""
mutable struct ExchangeProblem
    experiments::Vector{AbstractExperiment}
    model::AbstractModel
    integration::Union{Nothing,NamedTuple{(:peakppm, :noiseppm, :ppmwidth)}}
    savecalibration::Any
end
function ExchangeProblem(experiments, model, integration=nothing)
    return ExchangeProblem(experiments, model, integration, nothing)
end

"""
    FitResult

Result of a joint fit of exchange experiments. Fields are accessible via dot syntax.

# Fields
- `params`: fitted parameters with uncertainties (`ComponentArray{Measurement}`)
- `params_value`: fitted parameters as plain `Float64`, for simulation/plotting
- `params0`: initial parameters used for the fit
- `chi2`: chi-squared statistic
- `reduced_chi2`: chi-squared / degrees of freedom
- `cov`: parameter covariance matrix
- `nobs`: number of observations
- `nparams`: number of fitted (i.e. non-fixed) parameters
- `dof`: degrees of freedom
- `fixed`: flat indices of parameters held fixed during the fit (see `fit`)
- `prob`: the `ExchangeProblem` that was fitted
"""
struct FitResult
    params::ComponentArray
    params_value::ComponentArray
    params0::ComponentArray
    chi2::Float64
    reduced_chi2::Float64
    cov::Matrix{Float64}
    nobs::Int
    nparams::Int
    dof::Int
    fixed::Set{Int}
    prob::ExchangeProblem
end
