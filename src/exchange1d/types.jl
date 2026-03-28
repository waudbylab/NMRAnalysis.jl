abstract type AbstractModel end

abstract type AbstractExperiment end

struct ExchangeProblem
    experiments::Vector{AbstractExperiment}
    model::AbstractModel
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
- `nparams`: number of fitted parameters
- `dof`: degrees of freedom
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
    prob::ExchangeProblem
end