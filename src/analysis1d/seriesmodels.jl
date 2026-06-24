"""
    SeriesModel

Abstract supertype for the rule mapping a quantity series (vs an evolution parameter,
or between categorical slices) to derived parameters. Two concrete shapes:

- [`CurveFitModel`](@ref): a continuous fit-axis, fitted by nonlinear least squares.
- [`ContrastModel`](@ref): categorical slices combined algebraically (e.g. STD).

[`NoFitting`](@ref) simply carries the reduced quantities through (kinetics v1).
"""
abstract type SeriesModel end

"""
    CurveFitModel(func, param_names, estimate; xlabel, ylabel)

A parametric model fitted to `(x, y)` with `func(x, p)`. `estimate(x, y)` returns initial
parameters. Mirrors the `ParametricModel` pattern in GUI2D's `models.jl`.
"""
struct CurveFitModel <: SeriesModel
    func::Function
    param_names::Vector{String}
    estimate::Function
    xlabel::String
    ylabel::String
end

function CurveFitModel(func, param_names, estimate; xlabel="x", ylabel="Intensity")
    return CurveFitModel(func, collect(String, param_names), estimate, xlabel, ylabel)
end

"""
    NoFitting()

Pass the reduced quantities through without fitting a model (used for kinetics v1,
where the deliverable is intensity vs time).
"""
struct NoFitting <: SeriesModel end

# ---- concrete curve-fit models ------------------------------------------------

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

"""
    DampedSinusoidModel(; phase = :sine)

Damped sinusoid for nutation calibration: `A·sin(2π·ν·t)·exp(−R·t)` (or `cos` when
`phase = :cosine`), parameters `[A, ν, R]`. `ν` is the nutation frequency (Hz) from
which the 90° pulse length follows as `1/(4ν)`.
"""
function DampedSinusoidModel(; phase::Symbol=:sine)
    trig = phase === :cosine ? cos : sin
    function est(x, y)
        A = maximum(abs.(y))
        # rough frequency: assume ~half a period across the sampled range
        ν = 0.5 / (maximum(x) - minimum(x))
        return [A, ν, 1.0 / maximum(x)]
    end
    return CurveFitModel((x, p) -> @.(p[1] * trig(2π * p[2] * x) * exp(-p[3] * x)),
                         ["A", "ν", "R"],
                         est;
                         xlabel="Pulse duration / s")
end

# ---- contrast model (STD / waterLOGSY) ----------------------------------------

"""
    ContrastModel(; reference)

A categorical series model: within a group of slices it contrasts each non-reference
slice against the `reference` slice. The contrast value is `(I_ref − I_slice)/I_ref`
(the STD fraction). Used by the STD experiment; the 1D analogue of GUI2D's `hetnoe2d`.
"""
struct ContrastModel <: SeriesModel
    reference::Any
end

ContrastModel(; reference) = ContrastModel(reference)

# ---- fitting ------------------------------------------------------------------

"""
    fit_series(model, x, y) -> NamedTuple

Fit `model` to evolution values `x` and reduced quantities `y` (a
`Vector{Measurement}`). Returns `(; params, names, model, converged)` with `params` a
`Vector{Measurement}` (value ± standard error). The fit is noise-weighted using the
uncertainties carried by `y` (addressing the "scale fitting by noise" TODO).
"""
function fit_series(m::CurveFitModel, x::AbstractVector, y::AbstractVector)
    yv = Measurements.value.(y)
    yσ = Measurements.uncertainty.(y)
    p0 = m.estimate(x, yv)
    fit = if all(>(0), yσ)
        curve_fit(m.func, x, yv, 1.0 ./ yσ .^ 2, p0)
    else
        curve_fit(m.func, x, yv, p0)
    end
    params = coef(fit) .± stderror(fit)
    return (; params, names=m.param_names, model=m, converged=fit.converged)
end

function fit_series(::NoFitting, x::AbstractVector, y::AbstractVector)
    return (; params=Measurement{Float64}[], names=String[], model=NoFitting(),
            converged=true)
end
