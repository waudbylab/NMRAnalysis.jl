# Diffusion: Stejskal–Tanner analysis, giving D and (with solvent + temperature) rH.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    diffusion1d()
    diffusion1d(spec, gradients=nothing; coherence=SQ(H1), δ=nothing, Δ=nothing,
                σ=nothing, Gmax=nothing, regions=nothing, integration=nothing,
                prompt=isinteractive())

Analyse a diffusion experiment, fitting the Stejskal-Tanner equation to extract the
diffusion coefficient (and the hydrodynamic radius where the solvent and temperature are
known), then open the analysis window.

Called interactively, the acquisition parameters are read from the experiment and offered
for confirmation before the window opens, exactly as they must be: Bruker records the
gradient pulse length `δ` (`2·p30`), the diffusion delay `Δ` (`d20`) and the gradient
shape (`gpnam6`), but *not* the gradient ramp itself, so the initial and final gradient
strengths have to be supplied. Anything passed as an argument is taken as given and not
asked about.

# Arguments
- `gradients`: relative gradient strengths (0 to 1), one per spectrum. Supply these for
  any ramp other than a linear one, which is what the prompt constructs.
- `coherence`: the coherence whose gyromagnetic ratio applies (default `SQ(H1)`).
- `δ`: gradient pulse length, in seconds.
- `Δ`: diffusion delay, in seconds.
- `σ`: gradient shape factor (0.9 for the smoothed-square `SMSQ` shapes, 0.6366 for
  `SINE`).
- `Gmax`: maximum gradient strength, in T m⁻¹ (0.55 on a typical Bruker system).
- `regions`: integration regions to start from, instead of one on the tallest peak.
- `integration`: a `(; peakppm, noiseppm, ppmwidth)` triple, which skips the window and
  analyses that region directly.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.

# Examples
```julia
diffusion1d("106")                                    # asks for the gradient ramp
diffusion1d("106", 0.05:0.05:0.95)                    # ramp given explicitly
diffusion1d("106", g; δ=4e-3, Δ=0.1, σ=0.9, Gmax=0.55)  # nothing left to ask
```
"""
function diffusion1d(spec, gradients=nothing; coherence=SQ(H1), δ=nothing, Δ=nothing,
                     σ=nothing, Gmax=nothing, regions=nothing, integration=nothing,
                     prompt::Bool=isinteractive())
    spec = loadspec(spec)
    γ = gyromagneticratio(coherence)
    n = nplanesfromspec(spec)

    # Each parameter in turn: taken as given if supplied, otherwise read from the
    # acquisition parameters and offered for confirmation (the gradient list excepted -
    # Bruker does not record it, so there is nothing to confirm).
    prompt && println("Parsing experiment parameters...")
    p30 = acqusvalue(spec, :p, 30)
    δ = @something(δ,
                   1e-6 * ask("Gradient pulse length δ",
                              isnothing(p30) ? nothing : 2e6 * p30;
                              unit="µs", note=" (2*p30)", prompt))
    Δ = @something(Δ, ask("Diffusion delay Δ", acqusvalue(spec, :d, 20); unit="s",
                          note=" (d20)", prompt))
    gpnam = acqusvalue(spec, :gpnam, 6)
    σ = @something(σ, ask("Gradient shape factor σ", shapefactor(gpnam);
                          note=isnothing(gpnam) ? "" : " (gpnam6 = $gpnam)", prompt))
    Gmax = @something(Gmax, ask("Maximum gradient strength Gmax", 0.55; unit="T m⁻¹",
                                note=" (typical for Bruker systems)", prompt))
    gradients = @something(gradients, gradientramp(n; prompt))
    length(gradients) == n ||
        throw(ArgumentError("got $(length(gradients)) gradient strengths for $n spectra"))

    temp = acqusvalue(spec, :te)
    solvent = solventname(acqusvalue(spec, :solvent))

    ds = datasetfromspec(spec, [(; gradient=Float64(g)) for g in gradients])
    kw = (; γ, δ, Δ, σ, Gmax, temp, solvent)
    expt = isnothing(regions) ? DiffusionExperiment(ds; kw...) :
           DiffusionExperiment(ds; kw..., regions)
    return run1d(expt; integration)
end

diffusion1d(; kwargs...) = diffusion1d(askpath("diffusion experiment"); kwargs...)

"""
    gradientramp(n; prompt=true) -> Vector{Float64}

Ask for the gradient ramp, which Bruker does not store: the initial and final gradient
strengths (%) of a linear ramp over `n` spectra, matching TopSpin's `dosy` dialogue. Its
quadratic and exponential ramps are not generated here; pass those gradient strengths to
[`diffusion1d`](@ref) directly.
"""
function gradientramp(n::Integer; prompt::Bool=true)
    prompt ||
        throw(ArgumentError("the gradient ramp is not stored in the Bruker acquisition " *
                            "parameters - please pass `gradients` (relative strengths " *
                            "from 0 to 1, one per spectrum)"))
    println("The gradient ramp is not stored by Bruker, so it has to be entered here.")
    println("A linear ramp over $n spectra is assumed; for any other shape, " *
            "pass `gradients` instead.")
    g1 = ask("initial gradient strength"; unit="%", prompt)
    g2 = ask("final gradient strength"; unit="%", prompt)
    return collect(LinRange(0.01g1, 0.01g2, n))
end

"""Gradient shape factor from the Bruker shape name (as in the legacy routine), falling
back to 1.0 where the shape is unrecognised or was not recorded at all."""
function shapefactor(gpnam)
    (isnothing(gpnam) || length(gpnam) < 4) && return 1.0
    gpnam[1:4] == "SMSQ" && return 0.9
    gpnam[1:4] == "SINE" && return 0.6366
    return 1.0
end

"""The solvent as `viscosity` names it, or `nothing` for one it does not know (which
leaves the diffusion coefficient reported but not the hydrodynamic radius)."""
solventname(s) = s == "D2O" ? :d2o : (s == "H2O+D2O" ? :h2o : nothing)

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

# Own display names and units, not the shared PARAM_LABELS/PARAM_UNITS tables -
# everything about this experiment's presentation lives here. :D, :rH and :viscosity
# only ever appear in this file (StejskalTannerModel and postfit!, above); :rH has no
# entry in the labels table below (it stays the bare "rH" it already was) but does need
# its unit, which was previously sitting in the shared table for no reason but this.
const DIFFUSION_PARAM_LABELS = Dict(:D => "Diffusion coefficient",
                                    :viscosity => "η")
const DIFFUSION_PARAM_UNITS = Dict(:D => " ×10⁻¹⁰ m² s⁻¹",
                                   :rH => " Å",
                                   :viscosity => " mPa s")

function paramlabel(::DiffusionExperiment, name::Symbol)
    return get(DIFFUSION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::DiffusionExperiment, name::Symbol)
    return get(DIFFUSION_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
