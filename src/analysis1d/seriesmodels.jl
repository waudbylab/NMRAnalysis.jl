# Series models: the rule mapping a quantity series to derived parameters.
#
# Models used by two or more experiments live here; a model used by exactly one lives in
# that experiment's own `expt-*.jl` (StejskalTanner in expt-diffusion.jl, DampedSinusoid
# in expt-nutation.jl, Buildup and Contrast in expt-std.jl).

"""
    SeriesModel

Abstract supertype for the rule mapping a quantity series to derived parameters.

- [`CurveFitModel`](@ref): a continuous fit-axis, fitted by nonlinear least squares.
- [`NoFitting`](@ref): carry the reduced quantities through unfitted, for an experiment
  whose deliverable is the series itself (kinetics), or one that fits in `postfit!`
  because it must transform the series first (STD).
"""
abstract type SeriesModel end

"""
    CurveFitModel(func, paramnames, estimate; xlabel, ylabel)

A parametric model fitted to `(x, y)` with `func(x, p)`. `estimate(x, y)` returns initial
parameters. Mirrors the `ParametricModel` pattern in GUI2D's `models.jl`.
"""
struct CurveFitModel <: SeriesModel
    func::Function
    paramnames::Vector{String}
    estimate::Function
    xlabel::String
    ylabel::String
end

function CurveFitModel(func, paramnames, estimate; xlabel="x", ylabel="Intensity")
    return CurveFitModel(func, collect(String, paramnames), estimate, xlabel, ylabel)
end

"""
    NoFitting()

Pass the reduced quantities through without fitting a model (used for kinetics v1,
where the deliverable is intensity vs time).
"""
struct NoFitting <: SeriesModel end

# ---- shared curve-fit models ---------------------------------------------------

"""
    ExponentialModel()

Single exponential decay `A·exp(−R·t)`, parameters `[A, R]`. Suitable for T2/T1ρ
relaxation and TROSY/anti-TROSY decays.
"""
function ExponentialModel()
    return CurveFitModel((x, p) -> @.(p[1] * exp(-p[2] * x)),
                         ["A", "R"],
                         (x, y) -> [maximum(abs.(y)), 3.0 / maximum(x)];
                         xlabel="Time / s")
end

"""
    RecoveryModel()

Inversion/saturation recovery `A·(1 − C·exp(−R·t))`, parameters `[A, C, R]`.
"""
function RecoveryModel()
    return CurveFitModel((x, p) -> @.(p[1] * (1 - p[2] * exp(-p[3] * x))),
                         ["A", "C", "R"],
                         (x, y) -> [maximum(y), 2.0, 3.0 / maximum(x)];
                         xlabel="Time / s")
end

# ---- fitting ------------------------------------------------------------------

"""
    fitseries(model, x, y) -> NamedTuple

Fit `model` to evolution values `x` and reduced quantities `y` (a
`Vector{Measurement}`). Returns `(; params, names, model, converged)` with `params` a
`Vector{Measurement}` (value ± standard error). The fit is noise-weighted using the
uncertainties carried by `y` (addressing the "scale fitting by noise" TODO).
"""
function fitseries(m::CurveFitModel, x::AbstractVector, y::AbstractVector)
    yv = Measurements.value.(y)
    yσ = Measurements.uncertainty.(y)
    p0 = m.estimate(x, yv)
    fit = if all(>(0), yσ)
        curve_fit(m.func, x, yv, 1.0 ./ yσ .^ 2, p0)
    else
        curve_fit(m.func, x, yv, p0)
    end
    params = coef(fit) .± stderror(fit)
    return (; params, names=m.paramnames, model=m, converged=fit.converged)
end

function fitseries(::NoFitting, x::AbstractVector, y::AbstractVector)
    return (; params=Measurement{Float64}[], names=String[], model=NoFitting(),
            converged=true)
end
