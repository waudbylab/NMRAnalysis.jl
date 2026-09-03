# Diffusion: Stejskal–Tanner analysis, giving D and (with solvent + temperature) rH.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    diffusion1d(spec, gradients; coherence=SQ(H1), δ=nothing, Δ=nothing, σ=nothing,
                Gmax=0.55, regions=nothing, integration=nothing)

Analyse a diffusion experiment, fitting the Stejskal–Tanner equation to extract the
diffusion coefficient (and the hydrodynamic radius where the solvent and temperature are
known).

`gradients` is the vector of relative gradient strengths (0–1), one per plane. The
gradient pulse length `δ` (`2·p30`), diffusion delay `Δ` (`d20`) and shape factor `σ`
(from `gpnam6`) default to the values in the acquisition parameters.
"""
function diffusion1d(spec, gradients; coherence=SQ(H1), δ=nothing, Δ=nothing, σ=nothing,
                     Gmax=0.55, regions=nothing, integration=nothing)
    spec = _spec(spec)
    γ = gyromagneticratio(coherence)
    δ = isnothing(δ) ? acqus(spec, :p, 30) * 2.0 : δ
    Δ = isnothing(Δ) ? acqus(spec, :d, 20) : Δ
    σ = isnothing(σ) ? _shapefactor(acqus(spec, :gpnam, 6)) : σ

    temp = acqus(spec, :te)
    solvent = _solvent(acqus(spec, :solvent))

    ds = dataset_from_spec(spec, [(; gradient=Float64(g)) for g in gradients])
    kw = (; γ, δ, Δ, σ, Gmax, temp, solvent)
    expt = isnothing(regions) ? DiffusionExperiment(ds; kw...) :
           DiffusionExperiment(ds; kw..., regions)
    return run1d(expt; integration)
end

"""Gradient shape factor from the Bruker shape name (as in the legacy routine)."""
function _shapefactor(gpnam)
    length(gpnam) ≥ 4 || return 1.0
    gpnam[1:4] == "SMSQ" && return 0.9
    gpnam[1:4] == "SINE" && return 0.6366
    return 1.0
end

_solvent(s) = s == "D2O" ? :d2o : (s == "H2O+D2O" ? :h2o : nothing)

# ---- 2. type ------------------------------------------------------------------

"""
    DiffusionExperiment(dataset; γ, δ, Δ, σ, Gmax, temp=nothing, solvent=nothing, regions=…)

Fit integrated intensity vs relative gradient strength (`vars.gradient`, 0–1) to the
Stejskal–Tanner equation. The fitted `D` is the diffusion coefficient; `postfit!` adds the
hydrodynamic radius via the Stokes–Einstein relation when the solvent and temperature are
known.
"""
struct DiffusionExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
    temp::Union{Nothing,Float64}
    solvent::Any
end

function DiffusionExperiment(dataset::Dataset1D; γ, δ, Δ, σ, Gmax, temp=nothing,
                             solvent=nothing, regions=[defaultregion(dataset)])
    return DiffusionExperiment(dataset, collect(Region, regions),
                               StejskalTannerModel(; γ, δ, Δ, σ, Gmax),
                               isnothing(temp) ? nothing : Float64(temp), solvent)
end

# ---- 3. interface -------------------------------------------------------------

fitaxis(::DiffusionExperiment) = :gradient
primaryparam(::DiffusionExperiment) = :D

# ---- 4. science ---------------------------------------------------------------

"""
    StejskalTannerModel(; γ, δ, Δ, σ, Gmax)

Stejskal–Tanner diffusion decay `A·exp(−(γ·δ·σ·g·Gmax)²·(Δ−δ/3)·D)`, fitted against the
relative gradient strength `g` (0–1). Parameters are `[A, D]` with `D` in units of
10⁻¹⁰ m² s⁻¹ (as in the legacy routine), so the physical coefficient is `D · 1e-10`.
"""
function StejskalTannerModel(; γ, δ, Δ, σ, Gmax)
    k = (γ * δ * σ * Gmax)^2 * (Δ - δ / 3) * 1e-10
    return CurveFitModel((x, p) -> @.(p[1] * exp(-k * x^2 * p[2])),
                         ["A", "D"],
                         (x, y) -> [maximum(abs.(y)), 1.0];
                         xlabel="Relative gradient strength")
end

# D itself needs no post-fit: the fitted parameter is already the diffusion coefficient in
# the ×10⁻¹⁰ m² s⁻¹ that such coefficients are conventionally quoted in (and that the fit
# is conditioned on), so it is reported straight from `parameters`. Only the
# Stokes–Einstein radius is derived, and only when the solvent and temperature are known.
function postfit!(r::RegionResult, e::DiffusionExperiment)
    (isnothing(e.solvent) || isnothing(e.temp)) && return nothing
    D = param(r, :D) * 1e-10                       # m² s⁻¹
    η = viscosity(e.solvent, e.temp)               # mPa s
    kB = 1.38e-23
    setpost!(r, :viscosity, η)
    setpost!(r, :rH, kB * e.temp / (6π * η * 0.001 * D) * 1e10)   # Å
    return nothing
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::DiffusionExperiment) = "Diffusion"

resultlabels(::DiffusionExperiment) = ("Relative gradient strength",
                                        "Integrated intensity (a.u.)")

function spectruminfo(::DiffusionExperiment, vars::NamedTuple)
    return "gradient = $(round(100 * vars.gradient; digits=1))%"
end
