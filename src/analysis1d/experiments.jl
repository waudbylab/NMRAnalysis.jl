"""
    Experiment1D

Abstract supertype for 1D analyses. A concrete experiment is a thin composition that
supplies a dataset, a list of regions, a reduction, a series model, and the
fit-axis / grouping designation. The generic [`analyse`](@ref) pipeline does the rest;
experiments override [`postprocess`](@ref) to derive a global result (e.g. TRACT τc),
and may override `analyse` entirely for non-curve-fit shapes (e.g. STD contrast).

# Adding an experiment

One file per experiment, `expt-<name>.jl`, included at the foot of this file, laid out in
a fixed five-section order (see `expt-tract.jl` for the fullest example):

1. **entry point** — the exported `<name>1d(...)` function: load the spectra, pull the
   acquisition parameters, build the experiment, `run1d(expt; integration)`.
2. **type** — the `struct <Name>Experiment <: Experiment1D` and its keyword constructor.
3. **interface** — only the hooks that differ from the defaults below.
4. **science** — the series model, the derived-quantity maths, physics constants.
5. **presentation** — `windowtitle`, `result_labels`, `spectruminfo`, and friends.

Interface (with defaults):
- `dataset(e)`        — the `Dataset1D` (default: `e.dataset`)
- `regions(e)`        — `Vector{Region}`, length ≥ 1 (default: `e.regions`)
- `reduction(e)`      — `Reduction` (default `Integrate()`)
- `seriesmodel(e)`    — `SeriesModel` (default: `e.model`)
- `fitaxis(e)`        — `Symbol` naming the evolution variable (required)
- `groupcols(e)`      — `Tuple` of grouping variables (default `()`)
- `postprocess(e, results)` — global derived result (default `nothing`)
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
    series_results(e, [dataset, regions]; isfitting=true) -> Vector{SeriesResult}

Run the reduction and per-series curve fit for every region and grouping key. This is
the curve-fit pipeline shared by relaxation, TRACT, nutation and kinetics.

The `dataset`/`regions` arguments default to the experiment's own, but can be supplied
explicitly so the GUI can refit live against interactively-positioned regions and noise.
`isfitting=false` (the GUI's Fitting toggle switched off) substitutes [`NoFitting`](@ref)
for the experiment's own model, so the reduced quantities (`x`/`y`) still come through
for the plotted points, but no `curve_fit` call runs and `params`/`names` come back
empty - not merely a display toggle, an actual "don't fit" switch.
"""
series_results(e::Experiment1D) = series_results(e, dataset(e), regions(e))

function series_results(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    red = reduction(e)
    model = isfitting ? seriesmodel(e) : NoFitting()
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
    analyse(e, [dataset, regions]; isfitting=true) -> NamedTuple

Run the full analysis: `(; series, summary)` where `series` is a `Vector{SeriesResult}`
and `summary` is the experiment-specific global result (or `nothing`).

With `isfitting=false`, `postprocess` runs against an empty series list rather than
being skipped outright: several `postprocess` overrides (e.g. TRACT's) return a
`Vector`, and skipping straight to `nothing` would change `summary`'s type between
"fitting on" and "fitting off" - fatal for the GUI's `state[:result]` Observable, whose
element type is fixed by its first value. An empty series list keeps the same
(now-empty) `Vector` type either way.
"""
analyse(e::Experiment1D) = analyse(e, dataset(e), regions(e))

function analyse(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    series = series_results(e, ds, regs; isfitting)
    summary = postprocess(e, isfitting ? series : SeriesResult[])
    return (; series, summary)
end

# =============================================================================
# shared helpers
# =============================================================================

"""
    Integration(peakppm, noiseppm, ppmwidth)

The integration triple shared with `Exchange1D` (`prob.integration`): a peak position, a
noise position, and a common width, all in ppm (the noise region always has the same
width as the signal region — see [`reduce_region`](@ref)). Passing one to a top-level
entry point skips the GUI and analyses directly, so a previously-chosen region can be
replayed reproducibly from a script.
"""
const Integration = NamedTuple{(:peakppm, :noiseppm, :ppmwidth)}

regions_from(i) = [Region("signal", i.peakppm - i.ppmwidth / 2, i.peakppm + i.ppmwidth / 2)]

"""
    run1d(expt; integration=nothing)

Launch the GUI for `expt`, or - when an `integration` triple is supplied - skip the GUI
and return the analysis for that region directly.
"""
function run1d(expt::Experiment1D; integration=nothing)
    isnothing(integration) && return gui!(expt)
    ds = Dataset1D(dataset(expt).planes, Float64(integration.noiseppm))
    return analyse(expt, ds, regions_from(integration))
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
include("expt-std.jl")
include("expt-kinetics.jl")
