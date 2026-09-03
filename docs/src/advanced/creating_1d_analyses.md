# Creating New 1D Analyses

Every 1D analysis lives in a single self-contained file, `src/analysis1d/expt-<name>.jl`.
That file holds everything specific to the analysis — how the data is loaded, what is
fitted, what is derived from the fit, and how it is presented — so there is exactly one
place to go to change the science of an experiment, and one file to copy to add a new one.

The framework is deliberately the same shape as [`GUI2D`](creating_2d_analyses.md): an
abstract experiment type with defaulted hooks, a uniform result container, and a
visualisation strategy orthogonal to the experiment hierarchy. If you know one, you know
the other.

## The core idea

Every 1D-derived analysis is a collection of 1D traces, each tagged with the values of one
or more arrayed variables, reduced to quantities over named regions, then fitted against an
evolution parameter. An experiment is a thin composition over that shared machinery:

```
Dataset1D ──reduce──▶ quantity per plane ──group──▶ series ──fit──▶ parameters
                                                                       │
                                                              postfit! │ postfitglobal!
                                                                       ▼
                                                                 postparameters
```

## Data model

| type | meaning |
|---|---|
| `Trace(δ, y)` | one 1D spectrum: chemical shifts and intensities, plain vectors |
| `Planes(traces, vars)` | one row per spectrum; `vars[i]` is a `NamedTuple` of that spectrum's arrayed variables |
| `Region(label, lo, hi)` | a named ppm interval; zero width means a peak height |
| `Dataset1D(planes, noisecentre, label)` | the planes, the noise position, and where the data came from |

The analysis layer works on plain vectors and never sees an `NMRData`. The adapters in
`nmrdata.jl` (`loadspec`, `tracesfromspec`, `datasetfromspec`) are the only place NMRData
is touched.

## The file template

Each `expt-*.jl` is laid out in a fixed five-section order. `expt-tract.jl` exercises all
five and is the best worked example.

### 1. Entry point

The exported user-facing function. It loads the spectra, reads whatever it needs from the
acquisition parameters, builds the experiment, and hands off to `run1d`:

```julia
function myanalysis1d(spec; regions=nothing, integration=nothing)
    spec = loadspec(spec)
    times = something(annotation(spec, :myanalysis, :duration), acqus(spec, :vdlist))
    ds = datasetfromspec(spec, [(; time=Float64(t)) for t in times])
    expt = isnothing(regions) ? MyExperiment(ds) : MyExperiment(ds; regions)
    return run1d(expt; integration)
end
```

`run1d` opens the GUI, or — when an `integration = (; peakppm, noiseppm, ppmwidth)` triple
is supplied — skips it and analyses that region directly, so a chosen region can be
replayed from a script.

Any physics needed to interpret the acquisition parameters belongs here, in this file,
beside the science that uses it — not in a shared loader.

### 2. Type

```julia
struct MyExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function MyExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)])
    return MyExperiment(dataset, collect(Region, regions), ExponentialModel())
end
```

Fields named `dataset`, `regions` and `model` are picked up by the default accessors, so
no `dataset(e)`/`regions(e)`/`seriesmodel(e)` methods are needed. Override one only if your
experiment computes it rather than storing it.

### 3. Interface

Only what differs from the defaults:

| hook | default | purpose |
|---|---|---|
| `dataset(e)` | `e.dataset` | the data |
| `regions(e)` | `e.regions` | initial integration regions |
| `reduction(e)` | `Integrate()` | region × planes → quantity series |
| `seriesmodel(e)` | `e.model` | what is fitted |
| `fitaxis(e)` | *required* | the arrayed variable forming the x-axis |
| `groupcols(e)` | `()` | variables that split planes into separate series |
| `postfit!(r, e)` | none | quantities derived from one series |
| `postfitglobal!(results, e)` | none | quantities spanning several series |
| `primaryparam(e)` | `:A` | the headline quantity |
| `visualisationtype(e)` | `SeriesVisualisation()` | how results are drawn |

A *series* is the set of planes sharing all grouping variables, differing only in the
fit-axis. TRACT groups by `:which` (TROSY vs anti-TROSY), STD by `:sat`.

### 4. Science

Any model used by only this experiment is defined here; models shared by two or more
(`ExponentialModel`, `RecoveryModel`) live in `seriesmodels.jl`.

Derived quantities are computed in `postfit!` (one series) or `postfitglobal!` (several)
and recorded with `setpost!`:

```julia
function postfit!(r::RegionResult, ::MyExperiment)
    setpost!(r, :halflife, log(2) / param(r, :R))
    return nothing
end
```

`postfit!` is also where an experiment does its own fitting when the series must be
transformed first — STD contrasts each saturation series against the reference before
fitting the buildup curve, exactly as `cest2d` normalises against its reference plane. Use
`postfitglobal!` when the transformation needs another series, as STD's and TRACT's do.

Store each derived quantity **in the unit `PARAM_UNITS` names for it** (a 90° pulse in µs,
τc in ns), so one stored number serves both the summary and `results.csv`. Keys are unique
across all experiments — a symbol means one quantity in one unit everywhere — which is why
TRACT's cross-correlated rate is `:ηxy` and diffusion's viscosity is `:viscosity`.

### 5. Presentation

```julia
windowtitle(::MyExperiment) = "My analysis"
resultlabels(::MyExperiment) = ("Delay / s", "Integrated intensity (a.u.)")
spectruminfo(::MyExperiment, vars::NamedTuple) = "$(round(vars.time; digits=3)) s delay"
```

Optionally `resultxfactor` (display scaling for the x-axis) and `seriesnames` (a fixed
legend, where the groups are few and consistently named).

## Results

Every experiment fills the same container, `RegionResult` — the 1D analogue of GUI2D's
`Peak`:

- `parameters` — the fit's own output, keyed by `Symbol`
- `postparameters` — quantities derived from it by `postfit!`/`postfitglobal!`

Because the container is uniform, everything downstream is generic: the summary text, the
results panel and `results.csv` all read these two dictionaries without knowing which
analysis produced them. `analyse(expt)` returns a `Vector{RegionResult}`.

## Visualisation strategies

Presentation is a hierarchy *orthogonal* to `Experiment1D`, joined by
`visualisationtype(expt)` — the same arrangement as GUI2D's `VisualisationStrategy`. A
strategy implements three methods:

- `completeresultstate!(state, expt, ::V)` — build the Observables the panel reads
- `resultpanel!(gui, state, expt, ::V)` — the live GUI panel
- `plotresult!(ax, expt, result, label, i0, ::V)` — the static export

Live and export share one data getter, so what is saved is what was on screen. A strategy's
getter may return whatever shape suits it: only the three methods paired with it consume
it.

`SeriesVisualisation` (points, error bars, fitted line, one colour per group) is the
default and covers every experiment currently implemented. Declare a new strategy in your
own `expt-*.jl` only if that shape genuinely will not do. A strategy can also be held as a
*field* on the experiment, as GUI2D's `IntensityExperiment` does, if one experiment type
needs to be presented two different ways.

## Checklist for a new experiment

1. Write `src/analysis1d/expt-<name>.jl` following the five sections.
2. `include` it at the foot of `experiments.jl`.
3. Export the entry point and the experiment type in `Analysis1D.jl`, and re-export the
   entry point from `src/NMRAnalysis.jl`.
4. Add any new parameter symbols to `PARAM_UNITS` / `PARAM_LABELS` in `visualisation.jl`.
5. If it can run from a filename alone, register it in `Analysis1D.__init__` so `analyse`
   dispatches to it.

## Correspondence with GUI2D

| GUI2D | Analysis1D |
|---|---|
| `Experiment` | `Experiment1D` |
| `Peak` | `Region` (input) + `RegionResult` (output) |
| `SpecData` | `Planes` + `Dataset1D` |
| `peak.parameters` / `peak.postparameters` | `RegionResult.parameters` / `.postparameters` |
| `postfit!(peak, expt)` / `postfitglobal!(expt)` | `postfit!(r, expt)` / `postfitglobal!(results, expt)` |
| `primaryparam(expt)` | `primaryparam(expt)` |
| `preparespecdata` | each experiment's entry point + `nmrdata.jl` |
| `slicelabel(expt, idx)` | `spectruminfo(expt, vars)` |
| `peakinfotext(expt, idx)` | generic `resultsheader` / `secondarytext` |
| `experimentinfo(expt)` | `experimentinfo(expt)` (generic) |
| `VisualisationStrategy` | `ResultVisualisation` |
| `visualisationtype` / `completestate!` / `makepeakplot!` / `plot_peak!` | `visualisationtype` / `completeresultstate!` / `resultpanel!` / `plotresult!` |
| `files.jl` (`results.csv`, peak round-trip) | `files.jl` (`results.csv`, region round-trip) |
| `clustering.jl`, `mask!`, `simulate!` | — (no lineshape fitting in 1D) |
| `summaryplot` (parameter vs residue) | — (no residue axis in 1D) |

Deliberate differences: 1D has no `Parameter` type (there is no packed least-squares
vector, so a `Measurement` suffices), no fixed/moving split (regions never track between
spectra), and models are values rather than a type hierarchy — which means a 1D model
cannot carry its own `postfit!`, so derived science lives on the experiment.
