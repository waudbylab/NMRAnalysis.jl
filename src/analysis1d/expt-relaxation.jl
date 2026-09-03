# Relaxation: T1/T2 decay, or inversion/saturation recovery.
#
# Experiment files follow a fixed five-section order (see the module docs):
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    relaxation1d(spec; ir=nothing, tau=nothing, regions=nothing, integration=nothing)

Analyse a 1D relaxation experiment (T1/T2, or inversion recovery) from a pseudo-2D `spec`
(a path or an `NMRData`).

Relaxation delays come from `tau`, else the `relaxation.duration` annotation, else the
`vdlist`. The model defaults to the `relaxation.model` annotation
(`"inversion_recovery"` selects a recovery fit); pass `ir=true`/`false` to override.
"""
function relaxation1d(spec; ir=nothing, tau=nothing, regions=nothing, integration=nothing)
    spec = _spec(spec)
    times = something(tau, _ann(spec, :relaxation, :duration), acqus(spec, :vdlist))
    isrecovery = isnothing(ir) ?
                 (_ann(spec, :relaxation, :model) == "inversion_recovery") : ir
    ds = dataset_from_spec(spec, [(; time=Float64(t)) for t in times])
    expt = isnothing(regions) ? RelaxationExperiment(ds; ir=isrecovery) :
           RelaxationExperiment(ds; ir=isrecovery, regions)
    return run1d(expt; integration)
end

# ---- 2. type ------------------------------------------------------------------

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

function RelaxationExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)],
                              ir::Bool=false)
    return RelaxationExperiment(dataset, collect(Region, regions),
                                ir ? RecoveryModel() : ExponentialModel())
end

# ---- 3. interface -------------------------------------------------------------
# dataset/regions/seriesmodel come from the field-assuming defaults in experiments.jl

fitaxis(::RelaxationExperiment) = :time

# ---- 4. science ---------------------------------------------------------------
# The fitted rate is the deliverable; nothing further is derived from it.

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::RelaxationExperiment) = "Relaxation"

resultxfactor(e::RelaxationExperiment) = timescale(column(dataset(e).planes, :time))[1]

function result_labels(e::RelaxationExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end

spectruminfo(::RelaxationExperiment, vars::NamedTuple) = "$(round(vars.time; digits=3)) s delay"
