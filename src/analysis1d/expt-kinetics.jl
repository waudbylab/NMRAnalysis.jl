# Kinetics: intensity vs time over one or more named regions, possibly several runs.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    kinetics1d(spec, times; run=nothing, regions=nothing, model=NoFitting(), integration=nothing)

Track integrated intensity of one or more named regions against `times`, optionally tagged
by `run` for several time series.

`regions` defaults to a single peak-detected region; add further named regions of
interest interactively in the GUI (press `A`) rather than specifying them up front.
"""
function kinetics1d(spec, times::AbstractVector; run=nothing, regions=nothing,
                    model::SeriesModel=NoFitting(), integration=nothing)
    spec = _spec(spec)
    vars = if isnothing(run)
        [(; time=Float64(times[i])) for i in eachindex(times)]
    else
        [(; time=Float64(times[i]), run=run[i]) for i in eachindex(times)]
    end
    ds = dataset_from_spec(spec, vars)
    expt = isnothing(regions) ? KineticsExperiment(ds; model) :
           KineticsExperiment(ds; regions, model)
    return run1d(expt; integration)
end

# ---- 2. type ------------------------------------------------------------------

"""
    KineticsExperiment(dataset; regions=[defaultregion(dataset)], model=NoFitting())

Track integrated intensity of one or more named regions vs time (`vars.time`), grouped
by run (`vars.run`) when present. `regions` defaults to a single peak-detected region;
use the GUI (press `A`) to add further named regions of interest interactively rather
than specifying them up front. v1 defaults to `NoFitting` (the deliverable is the
intensity-vs-time trace); a kinetic `CurveFitModel` can be supplied to fit a model.
A future iteration will add an NMF reduction over selected ROIs.
"""
struct KineticsExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function KineticsExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)],
                            model::SeriesModel=NoFitting())
    return KineticsExperiment(dataset, collect(Region, regions), model)
end

# ---- 3. interface -------------------------------------------------------------

fitaxis(::KineticsExperiment) = :time
groupcols(e::KineticsExperiment) = hasvar(e.dataset.planes, :run) ? (:run,) : ()

# ---- 4. science ---------------------------------------------------------------
# v1 derives nothing beyond the intensity series itself.

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::KineticsExperiment) = "Kinetics"

result_labels(::KineticsExperiment) = ("Time", "Integrated intensity (a.u.)")

function spectruminfo(::KineticsExperiment, vars::NamedTuple)
    s = "$(round(vars.time; digits=3)) s"
    haskey(vars, :run) && (s *= " (run $(vars.run))")
    return s
end
