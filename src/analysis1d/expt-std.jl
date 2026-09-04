# STD: saturation-transfer difference, with buildup curves and an epitope map.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    std1d(spec, sat, tsat; regions=nothing, reference=:reference, excess=1.0, integration=nothing)

Analyse an STD experiment. `sat` and `tsat` are per-plane vectors giving the saturation
condition (one of which is the `reference`) and the saturation time. Reports STD fractions
per region, buildup curves where several saturation times are present, and an epitope map.

`regions` defaults to a single peak-detected region; add further named ligand regions
interactively in the GUI (press `A`) rather than specifying them up front.
"""
function std1d(spec, sat::AbstractVector, tsat::AbstractVector; regions=nothing,
               reference=:reference, excess::Real=1.0, integration=nothing)
    spec = loadspec(spec)
    vars = [(; sat=sat[i], tsat=Float64(tsat[i])) for i in eachindex(sat)]
    ds = datasetfromspec(spec, vars)
    expt = isnothing(regions) ? STDExperiment(ds; reference, excess) :
           STDExperiment(ds; regions, reference, excess)
    return run1d(expt; integration)
end

# ---- 2. type ------------------------------------------------------------------

"""
    STDExperiment(dataset; regions=[defaultregion(dataset)], reference=:reference, excess=1.0)

Saturation-transfer difference analysis. The planes carry two arrayed variables:

- `sat`  : the saturation condition — one designated `reference` (off-resonance) value
           plus one or more on-resonance saturation frequencies
           (e.g. `:reference, :methyl, :aromatic`).
- `tsat` : the saturation time (s).

For each ligand `region`, each non-reference saturation `sat`, and each `tsat`, the STD
fraction is contrasted against the reference at matching `tsat`:

    STD = (I_reference − I_sat) / I_reference

and the STD amplification factor `STD-AF = STD · excess`, where `excess = [L]/[P]`
(default 1, i.e. report STD%). When several `tsat` values are present for a
(region, sat) pair, a buildup curve `STD-AF(t) = STD-AF_max·(1 − exp(−k·t))` is fitted
and the initial slope `STD-AF₀ = STD-AF_max·k` reported (this removes the T1 bias that
makes raw STD% an unreliable epitope ranking). Finally an epitope map normalises STD-AF
across regions to the strongest signal.

This is a *contrast* experiment - the 1D analogue of GUI2D's `hetnoe2d`/`cest2d` - but it
runs through the ordinary pipeline: the reduction produces raw integrals per saturation
frequency, and [`postfitglobal!`](@ref) rewrites each non-reference series as its STD
fractions and fits the buildup curve, just as `postfit!(peak, ::CESTExperiment)` normalises
against its reference plane before fitting. It needs the *global* hook rather than the
per-series one only because in 1D the reference is a separate series (a distinct `sat`
value), where in 2D it is one of the peak's own planes.
"""
struct STDExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    reference::Any
    excess::Float64
end

function STDExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)],
                       reference=:reference, excess::Real=1.0)
    haskey(first(dataset.planes.vars), :sat) ||
        throw(ArgumentError("STD planes must carry a :sat variable"))
    haskey(first(dataset.planes.vars), :tsat) ||
        throw(ArgumentError("STD planes must carry a :tsat variable"))
    reference in column(dataset.planes, :sat) ||
        throw(ArgumentError("reference saturation $(repr(reference)) not present in :sat"))
    return STDExperiment(dataset, collect(Region, regions), reference, Float64(excess))
end

# ---- 3. interface -------------------------------------------------------------

# The reduction yields raw integrals per (region, sat); the contrast against the reference
# happens in `postfitglobal!`, which then fits the buildup curve, so no series model is
# fitted up front.
seriesmodel(::STDExperiment) = NoFitting()
fitaxis(::STDExperiment) = :tsat
groupcols(::STDExperiment) = (:sat,)
primaryparam(::STDExperiment) = :STD_AF

# ---- 4. science ---------------------------------------------------------------

"""
    BuildupModel()

STD buildup `STD-AF(t) = STD-AF_max·(1 − exp(−k·t))`, parameters `[STD_AF_max, k]`. The
reported epitope measure is the initial slope `STD-AF₀ = STD-AF_max·k`, which removes the
T1 bias that makes a raw STD% at a single saturation time an unreliable ranking.
"""
function BuildupModel()
    return CurveFitModel((x, p) -> @.(p[1] * (1 - exp(-p[2] * x))),
                         ["STD_AF_max", "k"],
                         (x, y) -> [maximum(y), 3.0 / maximum(x)];
                         xlabel="Saturation time / s", ylabel="STD-AF")
end

"""
    postfitglobal!(results, e::STDExperiment)

Turn the raw per-saturation integral series into STD series: contrast each non-reference
saturation against the reference at matching saturation time, fit the buildup curve where
there are enough times, then normalise across regions into an epitope map. The reference
series are dropped once they have been used - they carry no STD value of their own, and
the old readline routine did not report them either.
"""
function postfitglobal!(results::AbstractVector{RegionResult}, e::STDExperiment)
    isreference(r) = r.group.sat == e.reference
    for label in unique(r.region for r in results)
        rs = filter(r -> r.region == label, results)
        i = findfirst(isreference, rs)
        isnothing(i) && continue
        for r in rs
            isreference(r) || contrast!(r, rs[i], e.excess)
        end
    end
    filter!(!isreference, results)
    epitope!(results)
    return nothing
end

"""
    contrast!(r, reference, excess)

Rewrite `r` in place as the STD series for its saturation frequency: at each saturation
time present in both `r` and the `reference`, `STD = excess·(I_ref − I_sat)/I_ref`, with
replicates averaged on both sides. With three or more saturation times the buildup curve
is fitted, so that `r` ends up carrying an ordinary fit that the generic result panel and
summary render without knowing anything about STD.
"""
function contrast!(r::RegionResult, reference::RegionResult, excess::Real)
    at(s, τ) = mean(s.y[s.x .== τ])
    times = sort!(collect(Float64, intersect(unique(r.x), unique(reference.x))))
    stds = [excess * (at(reference, τ) - at(r, τ)) / at(reference, τ) for τ in times]
    r.x, r.y = times, stds
    if length(times) ≥ 3
        fit = fitseries(BuildupModel(), times, stds)
        r.parameters = OrderedDict{Symbol,Any}(Symbol(n) => p
                                               for (n, p) in zip(fit.names, fit.params))
        r.model = fit.model
        r.converged = fit.converged
        setpost!(r, :STD_AF0, param(r, :STD_AF_max) * param(r, :k))
    end
    # The epitope measure: the T1-unbiased initial slope where a buildup was fitted,
    # otherwise the STD fraction at the longest saturation time available.
    isempty(stds) ||
        setpost!(r, :STD_AF, get(r.postparameters, :STD_AF0, last(stds)))
    return r
end

"""
    epitope!(results)

Normalise the STD amplification factor across regions, separately for each saturation
frequency, to the strongest region - recording each as a percentage of it.
"""
function epitope!(results::AbstractVector{RegionResult})
    for sat in unique(r.group.sat for r in results)
        rs = filter(r -> r.group.sat == sat && haskey(r.postparameters, :STD_AF), results)
        isempty(rs) && continue
        values = [Measurements.value(r.postparameters[:STD_AF]) for r in rs]
        strongest = maximum(values)
        strongest == 0 && continue
        for (r, v) in zip(rs, values)
            setpost!(r, :relative, 100 * v / strongest)
        end
    end
    return results
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::STDExperiment) = "STD"

resultlabels(::STDExperiment) = ("Saturation time / s", "STD fraction")

spectruminfo(::STDExperiment, vars::NamedTuple) = "$(vars.sat), $(round(vars.tsat; digits=3)) s sat"

# Own display names and units, not the shared PARAM_LABELS/PARAM_UNITS tables -
# everything about this experiment's presentation lives here. :k, :STD_AF0, :STD_AF,
# :STD_AF_max and :relative only ever appear in this file (BuildupModel and
# postfitglobal!/epitope!, above).
const STD_PARAM_LABELS = Dict(:k => "Buildup rate",
                              :STD_AF0 => "STD-AF₀",
                              :STD_AF => "STD-AF",
                              :STD_AF_max => "STD-AF_max",
                              :relative => "epitope")
const STD_PARAM_UNITS = Dict(:k => " s⁻¹",
                             :relative => " %")

function paramlabel(::STDExperiment, name::Symbol)
    return get(STD_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::STDExperiment, name::Symbol)
    return get(STD_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
