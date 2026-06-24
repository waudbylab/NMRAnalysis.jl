"""
    Experiment1D

Abstract supertype for 1D analyses. A concrete experiment is a thin composition that
supplies a dataset, a list of regions, a reduction, a series model, and the
fit-axis / grouping designation. The generic [`analyse`](@ref) pipeline does the rest;
experiments override [`postprocess`](@ref) to derive a global result (e.g. TRACT τc),
and may override `analyse` entirely for non-curve-fit shapes (e.g. STD contrast).

Interface (with defaults):
- `dataset(e)`        — the `Dataset1D`
- `regions(e)`        — `Vector{Region}`, length ≥ 1
- `reduction(e)`      — `Reduction` (default `Integrate()`)
- `seriesmodel(e)`    — `SeriesModel`
- `fitaxis(e)`        — `Symbol` naming the evolution variable
- `groupcols(e)`      — `Tuple` of grouping variables (default `()`)
- `postprocess(e, results)` — global derived result (default `nothing`)
"""
abstract type Experiment1D end

reduction(::Experiment1D) = Integrate()
groupcols(::Experiment1D) = ()
postprocess(::Experiment1D, results) = nothing

"""
    SeriesResult

The fit of one series (one region × one grouping key).

# Fields
- `region`  : region label
- `group`   : `NamedTuple` of grouping values (empty if ungrouped)
- `x`       : evolution-parameter values (sorted)
- `y`       : reduced quantities (`Vector{Measurement}`)
- `params`  : fitted parameters (`Vector{Measurement}`)
- `names`   : parameter names
- `model`   : the series model used
- `converged` : whether the fit converged
"""
struct SeriesResult
    region::String
    group::NamedTuple
    x::Vector{Float64}
    y::Vector{Measurement{Float64}}
    params::Vector{Measurement{Float64}}
    names::Vector{String}
    model::Any
    converged::Bool
end

"""
    param(result, name) -> Measurement

Fitted value of parameter `name` from a `SeriesResult`.
"""
function param(r::SeriesResult, name::AbstractString)
    i = findfirst(==(name), r.names)
    isnothing(i) && throw(KeyError(name))
    return r.params[i]
end

"""
    series_results(e, [dataset, regions]) -> Vector{SeriesResult}

Run the reduction and per-series curve fit for every region and grouping key. This is
the curve-fit pipeline shared by relaxation, TRACT, nutation and kinetics.

The `dataset`/`regions` arguments default to the experiment's own, but can be supplied
explicitly so the GUI can refit live against interactively-positioned regions and noise.
"""
series_results(e::Experiment1D) = series_results(e, dataset(e), regions(e))

function series_results(e::Experiment1D, ds::Dataset1D, regs)
    red = reduction(e)
    model = seriesmodel(e)
    axis = fitaxis(e)
    results = SeriesResult[]
    for region in regs
        I = reduce_region(red, region, ds).I
        for (gkey, idx) in groupseries(ds.planes, groupcols(e))
            x = Float64[ds.planes.vars[i][axis] for i in idx]
            y = I[idx]
            perm = sortperm(x)
            x, y = x[perm], y[perm]
            fit = fit_series(model, x, y)
            push!(results,
                  SeriesResult(region.label, gkey, x, y, fit.params, fit.names, fit.model,
                               fit.converged))
        end
    end
    return results
end

"""
    analyse(e, [dataset, regions]) -> NamedTuple

Run the full analysis: `(; series, summary)` where `series` is a `Vector{SeriesResult}`
and `summary` is the experiment-specific global result (or `nothing`).
"""
analyse(e::Experiment1D) = analyse(e, dataset(e), regions(e))

function analyse(e::Experiment1D, ds::Dataset1D, regs)
    series = series_results(e, ds, regs)
    return (; series, summary=postprocess(e, series))
end

# =============================================================================
# Relaxation (T1/T2 decay, or inversion/saturation recovery)
# =============================================================================

"""
    RelaxationExperiment(dataset; ir=false)

Fit integrated intensity vs relaxation delay (`vars.time`) to a single exponential, or
to a recovery curve when `ir = true`. A single integration region (override via
`regions`).
"""
struct RelaxationExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function RelaxationExperiment(dataset::Dataset1D;
                              regions=[Region("signal",
                                              extrema_shift(dataset)...)],
                              ir::Bool=false)
    return RelaxationExperiment(dataset, collect(Region, regions),
                                ir ? RecoveryModel() : ExponentialModel())
end

dataset(e::RelaxationExperiment) = e.dataset
regions(e::RelaxationExperiment) = e.regions
seriesmodel(e::RelaxationExperiment) = e.model
fitaxis(::RelaxationExperiment) = :time

# =============================================================================
# TRACT (TROSY / anti-TROSY pair → τc)
# =============================================================================

"""
    TractExperiment(dataset; ωN, f, regions=…)

Fit TROSY and anti-TROSY decays (grouped by `vars.which ∈ {:trosy, :anti}`) and derive
the rotational correlation time τc from the cross-correlated relaxation-rate difference
ΔR = R(anti) − R(trosy). `ωN` is the ¹⁵N Larmor frequency (rad s⁻¹) and `f` the
dipole/CSA cross-correlation prefactor (see `tract_f`).
"""
struct TractExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    ωN::Float64
    f::Float64
end

function TractExperiment(dataset::Dataset1D; ωN, f,
                         regions=[Region("signal", extrema_shift(dataset)...)])
    return TractExperiment(dataset, collect(Region, regions), Float64(ωN), Float64(f))
end

dataset(e::TractExperiment) = e.dataset
regions(e::TractExperiment) = e.regions
seriesmodel(::TractExperiment) = ExponentialModel()
fitaxis(::TractExperiment) = :time
groupcols(::TractExperiment) = (:which,)

function postprocess(e::TractExperiment, results)
    summaries = NamedTuple[]
    for region in regions(e)
        rs = filter(r -> r.region == region.label, results)
        trosy = findfirst(r -> r.group.which == :trosy, rs)
        anti = findfirst(r -> r.group.which == :anti, rs)
        (isnothing(trosy) || isnothing(anti)) && continue
        Rtrosy = param(rs[trosy], "R")
        Ranti = param(rs[anti], "R")
        ηxy = (Ranti - Rtrosy) / 2
        τc = tract_tauc(e.f, e.ωN, ηxy)
        push!(summaries, (; region=region.label, Rtrosy, Ranti, ηxy, τc))
    end
    return summaries
end

"""
    tract_f(; B0, θ=17π/180) -> Float64

Dipole/CSA cross-correlation prefactor used in the TRACT τc relation, given the static
field `B0` (T). Constants follow the standard ¹⁵N–¹H amide treatment.
"""
function tract_f(; B0, θ=17 * π / 180)
    μ0 = 4π * 1e-7
    γH = 2.6752218744e8
    γN = -2.7126180e7
    ħ = 6.62607015e-34 / 2π
    rNH = 1.02e-10
    ΔδN = 160e-6
    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = B0 * γN * ΔδN / (3 * sqrt(2))
    return p * c * (3cos(θ)^2 - 1)
end

"""
    tract_tauc(f, ωN, ηxy) -> Float64

Rotational correlation time τc (ns) from the cross-correlated cross-relaxation rate
`ηxy`, by the analytic inversion of `ηxy = f·(4/5·τc + 3/5·τc/(1+(ωN·τc)²))` used in the
existing `tract` routine.
"""
function tract_tauc(f, ωN, ηxy)
    x = sqrt(21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy^2 * ωN^8 + 625 * f^2 * ηxy^4 * ωN^10)
    y = cbrt(1800 * f^2 * ηxy * ωN^4 + 125 * ηxy^3 * ωN^6 + 24 * sqrt(3) * x)
    τc = (5 * ηxy) / (24 * f) -
         (336 * f^2 * ωN^2 - 25 * ηxy^2 * ωN^4) / (24 * f * ωN^2 * y) + y / (24 * f * ωN^2)
    return 1e9 * τc
end

# =============================================================================
# Nutation calibration (pulse-length calibration)
# =============================================================================

"""
    NutationExperiment(dataset; phase=:sine, regions=…)

Fit integrated intensity vs pulse duration (`vars.duration`) to a damped sinusoid and
derive the 90° pulse length (`1/(4ν)`) and B₁ inhomogeneity (`R/2πν`). A height
(zero-width region at the observed signal) is the natural reduction but a finite region
works identically.
"""
struct NutationExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function NutationExperiment(dataset::Dataset1D; phase::Symbol=:sine,
                            regions=[Region("signal", extrema_shift(dataset)...)])
    return NutationExperiment(dataset, collect(Region, regions),
                              DampedSinusoidModel(; phase))
end

dataset(e::NutationExperiment) = e.dataset
regions(e::NutationExperiment) = e.regions
seriesmodel(e::NutationExperiment) = e.model
fitaxis(::NutationExperiment) = :duration

function postprocess(::NutationExperiment, results)
    return map(results) do r
        ν = param(r, "ν")
        R = param(r, "R")
        pulse90 = 1 / (4ν)
        inhomogeneity = R / (2π * ν)
        return (; region=r.region, ν, pulse90, inhomogeneity)
    end
end

# =============================================================================
# Kinetics (intensity vs time, possibly over multiple runs)
# =============================================================================

"""
    KineticsExperiment(dataset; regions, model=NoFitting())

Track integrated intensity of one or more named regions vs time (`vars.time`), grouped
by run (`vars.run`) when present. v1 defaults to `NoFitting` (the deliverable is the
intensity-vs-time trace); a kinetic `CurveFitModel` can be supplied to fit a model.
A future iteration will add an NMF reduction over selected ROIs.
"""
struct KineticsExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function KineticsExperiment(dataset::Dataset1D; regions, model::SeriesModel=NoFitting())
    return KineticsExperiment(dataset, collect(Region, regions), model)
end

dataset(e::KineticsExperiment) = e.dataset
regions(e::KineticsExperiment) = e.regions
seriesmodel(e::KineticsExperiment) = e.model
fitaxis(::KineticsExperiment) = :time
groupcols(e::KineticsExperiment) = hasvar(e.dataset.planes, :run) ? (:run,) : ()

# =============================================================================
# helpers
# =============================================================================

"""
    extrema_shift(dataset) -> (lo, hi)

Full chemical-shift span of the first plane — the default integration region when none
is supplied (the GUI will let the user narrow it).
"""
function extrema_shift(dataset::Dataset1D)
    δ = first(dataset.planes.traces).δ
    return (minimum(δ), maximum(δ))
end
