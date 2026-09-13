"""
    Experiment1D

Abstract supertype for 1D analyses. A concrete experiment is a thin composition supplying a
dataset, a list of regions, a series model, and the fit-axis and grouping designation; the
generic [`analyse`](@ref) pipeline does the rest. Only [`fitaxis`](@ref) is required; every
other hook has a default, and the defaults below assume fields named `dataset`, `regions`
and `model`.

One file per experiment, `expt-<name>.jl`, included at the foot of this file. See
`docs/src/advanced/creating_1d_analyses.md` for the file layout and the full list of hooks,
and `expt-tract.jl` for the fullest example.
"""
abstract type Experiment1D end

# Field-assuming defaults: an experiment that names its fields otherwise overrides one.
dataset(e::Experiment1D) = e.dataset
regions(e::Experiment1D) = e.regions
seriesmodel(e::Experiment1D) = e.model

groupcols(::Experiment1D) = ()

"""
    fitaxis(expt) -> Symbol

The arrayed variable forming the x-axis of every series: a relaxation delay, a gradient
strength, a pulse duration. Required of every experiment, there being no sensible default.
"""
function fitaxis end

"""
    integrate(region, expt[, dataset]) -> Vector{Measurement{Float64}}

Stage 1: reduce `region` to one measured quantity per plane. The default integrates;
an experiment that measures its regions some other way overrides this method.

`dataset` is separate from the experiment's own so the GUI can measure against
interactively-positioned regions and noise.
"""
integrate(region::Region, e::Experiment1D, ds::Dataset1D=dataset(e)) = integrate(region, ds)

"""
    SeriesResult

One measured series of a region: the planes sharing a grouping key, the quantity measured
in each, and the curve fitted through them.

`x` is the sorted evolution parameter, `y` the measured quantities, `planes` the dataset
plane each point came from (so the series file can name its source spectrum), and
`coefficients` the model's fitted coefficients in [`paramnames`](@ref) order, named onto
the region by [`liftparameters!`](@ref).
"""
struct SeriesResult
    group::NamedTuple
    x::Vector{Float64}
    y::Vector{Measurement{Float64}}
    planes::Vector{Int}
    model::Any
    coefficients::Vector{Measurement{Float64}}
    converged::Bool
end

"""
    RegionResult

Everything an analysis reports for one region: its `series`, one [`SeriesResult`](@ref) per
grouping key, and one flat set of `parameters`. The 1D analogue of GUI2D's `Peak`.

A series' fitted parameters are named apart by [`seriesname`](@ref) (`:R_trosy`, `:R_anti`)
and sit alongside whatever [`postfit!`](@ref) derives from them (`:tauc`), so one region is
one row of results. Values are untyped: a viscosity looked up from temperature is a plain
number where a fitted rate is a `Measurement`.

Parameters are stored in the unit [`paramunit`](@ref) names for them (a 90° pulse in µs, τc
in ns), not in SI, so one stored number serves both the summary and any tabular export.
"""
mutable struct RegionResult
    region::String
    series::Vector{SeriesResult}
    parameters::OrderedDict{Symbol,Any}
    postfitted::Bool
end

function RegionResult(region, series)
    return RegionResult(region, collect(SeriesResult, series),
                        OrderedDict{Symbol,Any}(), false)
end

"""
    seriesname(name, group) -> Symbol

The name a series' fitted parameter takes among its region's parameters: bare when the
region has one series, suffixed with the grouping values when it has several, so TRACT's
two decay rates become `:R_trosy` and `:R_anti`. Parameter names carry no underscore of
their own, so [`baseparam`](@ref) can strip the suffix again.
"""
function seriesname(name, group::NamedTuple)
    isempty(group) && return Symbol(name)
    return Symbol(name, "_", join(values(group), "_"))
end

"""
    baseparam(name) -> Symbol

The quantity a parameter name refers to, with any series suffix removed: `:R_trosy` → `:R`.
Labels and units are looked up under this name.
"""
function baseparam(name::Symbol)
    s = string(name)
    i = findfirst('_', s)
    return isnothing(i) ? name : Symbol(s[1:(i - 1)])
end

"""
    liftparameters!(result)

Name every series' fitted coefficients onto the region, via [`seriesname`](@ref). Runs as
each region is measured, so [`postfit!`](@ref) has only to add what it derives.
"""
function liftparameters!(r::RegionResult)
    for s in r.series
        for (name, value) in zip(paramnames(s.model), s.coefficients)
            r.parameters[seriesname(name, s.group)] = value
        end
    end
    return r
end

"""
    param(result, name) -> value

Value of parameter `name` (a `Symbol` or a `String`) for a region.
"""
param(r::RegionResult, name::Symbol) = r.parameters[name]
param(r::RegionResult, name::AbstractString) = param(r, Symbol(name))

"""
    isconverged(result) -> Bool

Whether every one of a region's series fitted successfully.
"""
isconverged(r::RegionResult) = all(s.converged for s in r.series)

"""
    setpost!(result, name, value)

Record a derived quantity on the region and mark it post-fitted.
"""
function setpost!(r::RegionResult, name::Symbol, value)
    r.parameters[name] = value
    r.postfitted = true
    return value
end

"""
    postfit!(result, expt)

Stage 2: derive the quantities reported for one region from its fitted series, recording
them with [`setpost!`](@ref) - nutation's 90° pulse from ν, diffusion's rH from D, TRACT's
τc from both rates. `result` holds every series of the region, so anything combining
conditions belongs here rather than in [`postfitglobal!`](@ref). The default derives
nothing.
"""
postfit!(::RegionResult, ::Experiment1D) = nothing

"""
    postfitglobal!(results, expt)

Stage 3: fit or derive quantities spanning every region of the analysis, and record them on
the relevant results. Runs after every `postfit!`. Mirrors GUI2D's `postfitglobal!(expt)`.
"""
postfitglobal!(::AbstractVector{RegionResult}, ::Experiment1D) = nothing

"""
    primaryparam(expt) -> Symbol

The experiment's headline quantity: the parameter a reader wants first, listed first
among the derived columns of any tabular export. Defaults to the amplitude `A` that every
`CurveFitModel` here fits. Mirrors GUI2D's `primaryparam(expt)`.
"""
primaryparam(::Experiment1D) = :A

"""
    seriesresults(e, [dataset, regions]; isfitting=true) -> Vector{RegionResult}

Measure and fit every region's series, for every grouping key, without the post-fit
stage. This is the pipeline shared by every curve-fit experiment.

The `dataset`/`regions` arguments default to the experiment's own, but can be supplied
explicitly so the GUI can refit live against interactively-positioned regions and noise.
`isfitting=false` (the GUI's Fitting toggle switched off) substitutes [`NoFitting`](@ref)
for the experiment's own model: the measured quantities still come through for the plotted
points, but no `curve_fit` call runs and the region's `parameters` come back empty.
"""
seriesresults(e::Experiment1D) = seriesresults(e, dataset(e), regions(e))

function seriesresults(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    model = isfitting ? seriesmodel(e) : NoFitting()
    axis = fitaxis(e)
    results = RegionResult[]
    for region in regs
        I = integrate(region, e, ds)
        series = SeriesResult[]
        for (gkey, idx) in groupseries(ds.planes, groupcols(e))
            x = Float64[ds.planes.vars[i][axis] for i in idx]
            y = I[idx]
            perm = sortperm(x)
            x, y, planes = x[perm], y[perm], idx[perm]
            fit = fitseries(model, x, y)
            push!(series,
                  SeriesResult(gkey, x, y, planes, fit.model, fit.params, fit.converged))
        end
        push!(results, liftparameters!(RegionResult(region.label, series)))
    end
    return results
end

"""
    analyse(e, [dataset, regions]; isfitting=true) -> Vector{RegionResult}

Run the full analysis: measure, fit, then post-fit. The return type does not depend on the
experiment or on whether anything was fitted, which is what lets the GUI's `state[:result]`
Observable keep a stable element type while it is still empty.
"""
analyse(e::Experiment1D) = analyse(e, dataset(e), regions(e))

function analyse(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    results = seriesresults(e, ds, regs; isfitting)
    isfitting || return results
    for r in results
        postfit!(r, e)
    end
    postfitglobal!(results, e)
    return results
end

# =============================================================================
# shared helpers
# =============================================================================

"""
    Integration(peakppm, noiseppm, ppmwidth)

The integration triple shared with `Exchange1D` (`prob.integration`): a peak position, a
noise position, and a common width, all in ppm. The noise region always takes the same
width as the signal region, so its own uncertainty matches (see [`integrate`](@ref)).
Passing one to an entry point skips the GUI and analyses directly.
"""
const Integration = NamedTuple{(:peakppm, :noiseppm, :ppmwidth)}

regionsfrom(i) = [Region("signal", i.peakppm - i.ppmwidth / 2, i.peakppm + i.ppmwidth / 2)]

"""
    run1d(expt; integration=nothing, call=nothing) -> Vector{RegionResult}

Launch the GUI for `expt` and return the results standing when its window is closed, or -
when an `integration` triple is supplied - skip the GUI and return the analysis for that
region directly.

`call` is the [`AnalysisCall`](@ref) the entry point recorded, carried through to
`summary.txt`; it is unused on the scripted path, which writes no files. Both paths return
the same thing. [`gui!`](@ref) returns the whole GUI state, of which this is one entry.
"""
function run1d(expt::Experiment1D; integration=nothing, call=nothing)
    isnothing(integration) && return gui!(expt; call)[:result][]
    d = dataset(expt)
    ds = Dataset1D(d.planes, Float64(integration.noiseppm), d.label, d.sources)
    return analyse(expt, ds, regionsfrom(integration))
end

"""
    defaultregion(dataset; label="signal", width=defaultregionwidth(...)) -> Region

A sensible default single integration region: `width` ppm wide (2% of the spectral
width by default) centred on the tallest peak (by absolute intensity) in the first
plane. Used as the default `regions` for every experiment with a single signal; the GUI
lets the user reposition and resize it, or add further regions.
"""
function defaultregion(dataset::Dataset1D; label="signal",
                       width=defaultregionwidth(first(dataset.planes.traces).δ))
    t = first(dataset.planes.traces)
    peak = t.δ[argmax(abs.(t.y))]
    return Region(label, peak - width / 2, peak + width / 2)
end

# =============================================================================
# implementations - one file per experiment
# =============================================================================

include("expt-relaxation.jl")
include("expt-tract.jl")
include("expt-nutation.jl")
include("expt-diffusion.jl")
include("expt-kinetics.jl")
