"""
    Experiment1D

Abstract supertype for 1D analyses. A concrete experiment is a thin composition that
supplies a dataset, a list of regions, a reduction, a series model, and the
fit-axis / grouping designation. The generic [`analyse`](@ref) pipeline does the rest;
experiments derive further quantities in [`postfit!`](@ref) / [`postfitglobal!`](@ref),
and may override `analyse` entirely for non-curve-fit shapes.

# Adding an experiment

One file per experiment, `expt-<name>.jl`, included at the foot of this file, laid out in
a fixed five-section order (see `expt-tract.jl` for the fullest example):

1. **entry point** — the exported `<name>1d(...)` function: load the spectra, pull the
   acquisition parameters, build the experiment, `run1d(expt; integration)`.
2. **type** — the `struct <Name>Experiment <: Experiment1D` and its keyword constructor.
3. **interface** — only the hooks that differ from the defaults below.
4. **science** — the series model, the derived-quantity maths, physics constants.
5. **presentation** — `windowtitle`, `resultlabels`, `spectruminfo`, and friends.

Interface (with defaults):
- `dataset(e)`        — the `Dataset1D` (default: `e.dataset`)
- `regions(e)`        — `Vector{Region}`, length ≥ 1 (default: `e.regions`)
- `reduction(e)`      — `Reduction` (default `Integrate()`)
- `seriesmodel(e)`    — `SeriesModel` (default: `e.model`)
- `fitaxis(e)`        — `Symbol` naming the evolution variable (required)
- `groupcols(e)`      — `Tuple` of grouping variables (default `()`)
- `postfit!(r, e)`    — derived quantities from one series (default: none)
- `postfitglobal!(results, e)` — derived quantities spanning series (default: none)
- `primaryparam(e)`   — the headline quantity (default `:A`)
"""
abstract type Experiment1D end

# Field-assuming defaults, in the style of GUI2D's `nslices(expt)`. An experiment whose
# fields are named `dataset`/`regions`/`model` needs none of these three methods; one
# that computes them (or names them differently) overrides just the one it needs.
dataset(e::Experiment1D) = e.dataset
regions(e::Experiment1D) = e.regions
seriesmodel(e::Experiment1D) = e.model

reduction(::Experiment1D) = Integrate()
groupcols(::Experiment1D) = ()

"""
    RegionResult

The analysis of one series: one region × one grouping key. The 1D analogue of GUI2D's
`Peak`, and like it the *uniform* container every experiment fills, however different
their science.

# Fields
- `region`    : region label
- `group`     : `NamedTuple` of grouping values (empty if ungrouped)
- `x`         : evolution-parameter values (sorted)
- `y`         : reduced quantities (`Vector{Measurement}`)
- `parameters`     : the curve fit's own parameters (`A`, `R`, `ν`, …)
- `postparameters` : quantities derived from them by [`postfit!`](@ref) /
  [`postfitglobal!`](@ref) (`τc`, `pulse90`, `rH`, …)
- `model`     : the series model used
- `converged` : whether the fit converged
- `postfitted`: whether a derived result has been recorded

Both parameter dictionaries are `Symbol`-keyed, mirroring `peak.parameters` /
`peak.postparameters`, so generic consumers (the summary formatter, `primaryparam`) read
them by name without knowing the experiment. Values are untyped because not every derived
quantity carries an uncertainty - a solvent viscosity looked up from temperature is a
plain number, where a fitted rate is a `Measurement`.

Derived quantities are stored **in the unit named by `PARAM_UNITS`** (a 90° pulse in µs,
τc in ns), not in SI, so that one stored number serves both the summary and any later
tabular export.
"""
mutable struct RegionResult
    region::String
    group::NamedTuple
    x::Vector{Float64}
    y::Vector{Measurement{Float64}}
    parameters::OrderedDict{Symbol,Any}
    postparameters::OrderedDict{Symbol,Any}
    model::Any
    converged::Bool
    postfitted::Bool
end

"""
    param(result, name) -> value

Fitted value of parameter `name` (a `Symbol` or a `String`) from a `RegionResult`.
"""
param(r::RegionResult, name::Symbol) = r.parameters[name]
param(r::RegionResult, name::AbstractString) = param(r, Symbol(name))

"""
    setpost!(result, name, value)

Record a derived quantity on `result` and mark it post-fitted. The `postfit!` /
`postfitglobal!` counterpart of writing into `peak.postparameters` and setting
`peak.postfitted[]` in GUI2D.
"""
function setpost!(r::RegionResult, name::Symbol, value)
    r.postparameters[name] = value
    r.postfitted = true
    return value
end

"""
    postfit!(result, expt)

Derive further quantities from one series' own fitted parameters and record them with
[`setpost!`](@ref) - nutation's 90° pulse length from ν, diffusion's rH from D. The
default does nothing (relaxation and kinetics derive nothing). Mirrors GUI2D's
`postfit!(peak, expt)`.
"""
postfit!(::RegionResult, ::Experiment1D) = nothing

"""
    postfitglobal!(results, expt)

Derive quantities that need more than one series - TRACT's τc, which combines the TROSY
and anti-TROSY rates for a region - and record them on the relevant results. Runs after
every `postfit!`. Mirrors GUI2D's `postfitglobal!(expt)`.
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

Run the reduction and per-series curve fit for every region and grouping key, without the
post-fit stage. This is the pipeline shared by every curve-fit experiment.

The `dataset`/`regions` arguments default to the experiment's own, but can be supplied
explicitly so the GUI can refit live against interactively-positioned regions and noise.
`isfitting=false` (the GUI's Fitting toggle switched off) substitutes [`NoFitting`](@ref)
for the experiment's own model, so the reduced quantities (`x`/`y`) still come through
for the plotted points, but no `curve_fit` call runs and `parameters` comes back empty -
not merely a display toggle, an actual "don't fit" switch.
"""
seriesresults(e::Experiment1D) = seriesresults(e, dataset(e), regions(e))

function seriesresults(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    red = reduction(e)
    model = isfitting ? seriesmodel(e) : NoFitting()
    axis = fitaxis(e)
    results = RegionResult[]
    for region in regs
        I = reduceregion(red, region, ds).I
        for (gkey, idx) in groupseries(ds.planes, groupcols(e))
            x = Float64[ds.planes.vars[i][axis] for i in idx]
            y = I[idx]
            perm = sortperm(x)
            x, y = x[perm], y[perm]
            fit = fitseries(model, x, y)
            parameters = OrderedDict{Symbol,Any}(Symbol(n) => p
                                                 for (n, p) in zip(fit.names, fit.params))
            push!(results,
                  RegionResult(region.label, gkey, x, y, parameters,
                               OrderedDict{Symbol,Any}(), fit.model, fit.converged, false))
        end
    end
    return results
end

"""
    analyse(e, [dataset, regions]; isfitting=true) -> Vector{RegionResult}

Run the full analysis: reduce, fit, then post-fit. The return type does not depend on the
experiment or on whether anything was fitted, so the GUI's `state[:result]` Observable has
a stable element type even when it starts out empty - which the old
`(; series, summary)` shape did not, `summary` being `nothing` for some experiments and a
`Vector` for others.
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
noise position, and a common width, all in ppm (the noise region always has the same
width as the signal region — see [`reduceregion`](@ref)). Passing one to a top-level
entry point skips the GUI and analyses directly, so a previously-chosen region can be
replayed reproducibly from a script.
"""
const Integration = NamedTuple{(:peakppm, :noiseppm, :ppmwidth)}

regionsfrom(i) = [Region("signal", i.peakppm - i.ppmwidth / 2, i.peakppm + i.ppmwidth / 2)]

"""
    run1d(expt; integration=nothing)

Launch the GUI for `expt`, or - when an `integration` triple is supplied - skip the GUI
and return the analysis for that region directly.
"""
function run1d(expt::Experiment1D; integration=nothing)
    isnothing(integration) && return gui!(expt)
    d = dataset(expt)
    ds = Dataset1D(d.planes, Float64(integration.noiseppm), d.label)
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
