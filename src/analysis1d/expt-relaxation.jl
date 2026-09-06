# Relaxation: T1/T2 decay, or inversion/saturation recovery.
#
# Experiment files follow a fixed five-section order (see the module docs):
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    relaxation1d(spec; tau=nothing, model=nothing, regions=nothing, integration=nothing,
                 prompt=isinteractive())

Analyse a 1D relaxation experiment (T1/T2, or inversion recovery) from a pseudo-2D `spec`
(a path or an `NMRData`), then open the analysis window.

Relaxation delays are taken from `tau` if given, else from the `relaxation.duration`
annotation, else from the `vdlist`; if none of those is available you are asked for them.
The model is taken from `model` if given, else from the `relaxation.model` annotation, and
otherwise you are asked to choose one. Neither annotation is required.

# Arguments
- `tau`: relaxation delays, in seconds, one per spectrum.
- `model`: `:exponential` for a decay `A exp(-Rt)`, `:recovery` for an
  inversion/saturation recovery `A [1 - C exp(-Rt)]`, or any [`SeriesModel`](@ref).
- `regions`: integration regions to start from, instead of one on the tallest peak.
- `integration`: a `(; peakppm, noiseppm, ppmwidth)` triple, which skips the window and
  analyses that region directly.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.

# Examples
```julia
relaxation1d("11")                                   # annotations or vdlist, or asks
relaxation1d("11"; tau=[0.01, 0.05, 0.1, 0.5])       # delays given explicitly
relaxation1d("11"; model=:recovery)                  # inversion recovery
```
"""
function relaxation1d(spec; tau=nothing, model=nothing, regions=nothing,
                      integration=nothing, prompt::Bool=isinteractive())
    spec = loadspec(spec)
    times = @something(tau,
                       annotation(spec, :relaxation, :duration),
                       acqusvalue(spec, :vdlist),
                       askvector("relaxation delays", nplanesfromspec(spec); unit="s",
                                 prompt))
    model = @something(model,
                       relaxationmodel(annotation(spec, :relaxation, :model)),
                       askrelaxationmodel(; prompt))
    ds = datasetfromspec(spec, [(; time=Float64(t)) for t in times])
    expt = isnothing(regions) ? RelaxationExperiment(ds; model) :
           RelaxationExperiment(ds; model, regions)
    return run1d(expt; integration)
end

"""
    relaxationmodel(annotation) -> Symbol or nothing

The model named by a `relaxation.model` annotation, or `nothing` where there is no
annotation to read - so the entry point can fall through to asking.
"""
relaxationmodel(::Nothing) = nothing
relaxationmodel(s) = string(s) == "inversion_recovery" ? :recovery : :exponential

"""Offer the choice of relaxation model. Without prompting, a plain decay is assumed."""
function askrelaxationmodel(; prompt::Bool=true)
    choice = askchoice("Which relaxation model should be fitted?",
                       ["Exponential decay:  I(t) = A exp(-Rt)",
                        "Inversion recovery: I(t) = A [1 - C exp(-Rt)]"]; prompt)
    return choice == 2 ? :recovery : :exponential
end

# ---- 2. type ------------------------------------------------------------------

"""
    RelaxationExperiment(dataset; model=:exponential, regions=…)

Fit integrated intensity vs relaxation delay (`vars.time`) to a single exponential
(`model = :exponential`) or to a recovery curve (`model = :recovery`); any
[`SeriesModel`](@ref) may also be passed directly. A single integration region (override
via `regions`).
"""
struct RelaxationExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function RelaxationExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)],
                              model=:exponential)
    return RelaxationExperiment(dataset, collect(Region, regions),
                                relaxationseriesmodel(model))
end

"""
    relaxationseriesmodel(model) -> SeriesModel

The curve fitted for a named relaxation model. A `SeriesModel` passes straight through, so
an unusual decay can be supplied directly without a name having to be invented for it.
"""
relaxationseriesmodel(model::SeriesModel) = model
function relaxationseriesmodel(model::Symbol)
    model === :exponential && return ExponentialModel()
    model === :recovery && return RecoveryModel()
    throw(ArgumentError("unknown relaxation model :$model " *
                        "(expected :exponential or :recovery)"))
end

# ---- 3. interface -------------------------------------------------------------
# dataset/regions/seriesmodel come from the field-assuming defaults in experiments.jl

fitaxis(::RelaxationExperiment) = :time
primaryparam(::RelaxationExperiment) = :R

# ---- 4. science ---------------------------------------------------------------
# The fitted rate is the deliverable; nothing further is derived from it.

"""
    RecoveryModel()

Inversion/saturation recovery `A·(1 − C·exp(−R·t))`, parameters `[A, C, R]`. Used only
here: TRACT always fits a plain exponential (`ExponentialModel`, shared in
`seriesmodels.jl` since both experiments use it), never a recovery curve.
"""
function RecoveryModel()
    return CurveFitModel((x, p) -> @.(p[1] * (1 - p[2] * exp(-p[3] * x))),
                         ["A", "C", "R"],
                         (x, y) -> [maximum(y), 2.0, 3.0 / maximum(x)];
                         xlabel="Time / s")
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::RelaxationExperiment) = "Relaxation"

resultxfactor(e::RelaxationExperiment) = timescale(column(dataset(e).planes, :time))[1]

function resultlabels(e::RelaxationExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end

spectruminfo(::RelaxationExperiment, vars::NamedTuple) = "$(round(vars.time; digits=3)) s delay"

# Own display names, not the shared PARAM_LABELS table - everything about this
# experiment's presentation lives here, not scattered into a module-wide table. :R
# genuinely means "relaxation rate" for both T1/T2 and recovery fits (TRACT overrides
# it too, for the same reason - see its own note); :C only ever appears here, since
# RecoveryModel, above, is this experiment's alone.
const RELAXATION_PARAM_LABELS = Dict(:R => "Relaxation rate",
                                     :C => "Recovery factor")

function paramlabel(::RelaxationExperiment, name::Symbol)
    return get(RELAXATION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
