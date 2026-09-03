# 1D Analysis Framework — Design & Implementation Plan

This document captures the design agreed in discussion for a unified 1D analysis
framework, and tracks implementation status.

## Philosophy

Many lightweight tools for particular analyses. Make 90% of the most common work
easier to use; leave the 10% (full lineshape deconvolution, CORCEMA-ST, numerical
Bloch–McConnell, NUS) alone. Each experiment should be a *thin composition* over
shared machinery, not a bespoke script.

## The core idea

Every 1D-derived analysis is **a collection of 1D traces, each tagged with values
of one or more arrayed variables**, reduced to quantities over named regions, then
fitted against an evolution parameter (or contrasted between categorical slices).

The differences between experiments are confined to three composition slots plus an
optional global post-fit:

1. **Reduction** — region × planes → named quantity series. v1: `Integrate` (a height
   is just a zero-width region; nutation needs no special case). Future: NMF for kinetics.
2. **Series model** — quantity vs evolution parameter → derived parameters. Two shapes:
   - **curve-fit** (continuous axis): exponential, recovery, damped sinusoid, …
   - **contrast** (categorical slices): STD = (I_ref − I_sat)/I_ref, waterLOGSY sign.
     (1D analogue of the existing `hetnoe2d` reference/saturated experiment.)
3. **Visualisation strategy** — usually derived from the series model; overridable.

Plus:

- **Post-fit** (`postfit!` / `postfitglobal!`) — derive further quantities from a fit and
  record them on the result (TRACT τc from ΔR₂; STD epitope normalisation), the same hooks
  GUI2D uses. `postfitglobal!` is the one to use when the derivation spans several series.
- **Noise region and the region list are universal** and live in the shared base.

## Data model

- `Trace(δ, y)` — one 1D spectrum: chemical-shift axis + intensities. Plain vectors, so
  the analysis layer is GUI-independent (the agreed "keep the science pure" split).
- `Planes(traces, vars)` — long format: one row per spectrum, `vars[i]` a `NamedTuple`
  of the arrayed variables (`:time`, `:which`, `:sat`, `:tsat`, `:run`, …). The name is
  `planes`, echoing the 2D `plane`/`slice`/`nslices` vocabulary.
- `Region(label, lo, hi)` — a named ppm interval; `lo == hi` ⇒ height.
- `Dataset1D(planes, noise)` — the planes plus the universal noise region.

Series = group plane indices by all `vars` columns except the fit-axis. This generalises
R1ρ's `onresseries` to arbitrary coordinates; TRACT groups by `:which`, kinetics by `:run`.

## Reactivity / GUI substrate (later phases)

- Keep the **science pure** (functions over `Trace`/`Planes`), Observables/ComputeGraph
  only at the GUI edge. Cost decides the binding: cheap/live derivations are reactive;
  the expensive, judgement-bearing fit (exchange) is button-triggered.
- Adopt Makie's **ComputeGraph (compute pipeline)** as the substrate for the *new* 1D
  GUI — lazy, cached, consistently resolved — which removes the eager-recompute and
  `fit_generation`-counter pains seen in GUI2D. Prototype it on the Tier-2 exchange popup
  first (highest value, new code, contained blast radius); migrate Window 1 if it pays off.
- Two windows: Window 1 = spectral overlay + auto-sliders (one per `vars` column with >1
  value) + region list + live Tier-1 result panel. Window 2 (opt-in per experiment) =
  exchange/global model fit. Both are views onto the same nodes; only the Fit button runs
  the heavy solve.

## Relationship to GUI2D

1D is GUI2D with the entire lineshape-fitting layer removed (no clustering, masking,
simulation, or threaded cancellable fit — reduction is just integration) and a
coordinate/second-window layer added. Reused wholesale: the `Parameter` idea, the
`models.jl` parametric-model pattern, the `postfit!`/`postfitglobal!` hooks, the
`VisualisationStrategy` trait + composition, state-as-reactive-nodes, `summaryplot`.
Genuinely new: `planes` (multi-axis navigation) and the exchange second window.

Module layout deliberately mirrors GUI2D so a shared `AnalysisCore` can later be lifted
out of both (start parallel, refactor toward shared core once 1D is proven). As of the
restructure below this is no longer approximate: the two modules share the same file
layout, the same `postfit!`/`postfitglobal!`/`primaryparam` hooks, the same uniform
`parameters`/`postparameters` result container, and the same
visualisation-strategy composition. See `docs/src/advanced/creating_1d_analyses.md` for the
full correspondence table.

## Experiments (this iteration)

| experiment | regions | plane vars | fit-axis | group | reduction | series model | global |
|---|---|---|---|---|---|---|---|
| relaxation | 1 | `time` | `time` | – | Integrate | Exponential / Recovery | – |
| TRACT | 1 | `time, which` | `time` | `which` | Integrate | Exponential ×2 | τc(ΔR₂) inline |
| nutation | 1 | `duration` | `duration` | – | Integrate | DampedSinusoid | 90° pulse |
| diffusion | 1 | `gradient` | `gradient` | – | Integrate | Stejskal–Tanner | D, rH |
| STD | N named | `sat, tsat` | `tsat` | `sat` | Integrate | Contrast (+ buildup) | epitope |
| kinetics | N named | `time, run` | `time` | `run` | Integrate | NoFitting (v1) | – |

Entry points take the names of the routines they replace (`relaxation1d`, `tract`,
`calibration1d`, `diffusion1d`, plus `std1d`, `kinetics1d`). Each opens the GUI; passing
`integration = (; peakppm, noiseppm, ppmwidth)` skips it and analyses that region
directly, so a chosen region can be replayed from a script. The triple is deliberately
the same one `Exchange1D` stores as `prob.integration`.

### STD details

STD is the richest case and is designed for:
- minimal: reference + one saturation, single `tsat` → STD% per region.
- multiple saturation frequencies (`sat ∈ {reference, methyl, aromatic, …}`) → STD% per
  (region, sat), each non-reference category contrasted against the reference at matching `tsat`.
- buildup: multiple `tsat` → fit STD-AF(tsat) = STD-AF_max·(1 − exp(−k·tsat)) per
  (region, sat); report the **initial slope** STD-AF₀ = STD-AF_max·k (removes T1 bias).
- epitope: normalise STD% across regions to the strongest → relative %.

Leave alone (10%): CORCEMA-ST relaxation-matrix epitope quantification.

## Implementation status

Phase 1 — **analysis core + 5 experiments**  ✓
- [x] `Trace` / `Planes` / `Region` / `Dataset1D`
- [x] `Integrate` reduction with noise propagation (Measurements), height = zero width
- [x] series models: Exponential, Recovery, DampedSinusoid, NoFitting; Contrast (STD)
- [x] grouping + curve-fit pipeline (noise-weighted) → `RegionResult`
- [x] experiments: Relaxation, TRACT (τc), Nutation (90°), STD (multi-freq + buildup + epitope), Kinetics
- [x] entry points `relaxation1d`, `tract`, `calibration1d`, `diffusion1d`, `std1d`,
      `kinetics1d`, each opening the GUI or honouring an `integration` triple

Phase 2 — **interactive GUI (Window 1)**  ✓
- [x] spectral overlay (all planes + current), draggable region(s) + noise region
- [x] plane slider, region selector, width boxes, fitting toggle, save
- [x] live refit via `lift`; result panel dispatched per experiment (`result_plotdata`)
- [x] launched with `gui!(expt)`; CairoMakie export of fit + summary
- [ ] ComputeGraph substrate (deferred; prototype on the exchange popup first)
- [ ] auto-generated per-axis sliders for multi-coordinate datasets

Phase 3 — **replace the readline routines + hook into dispatch**  ✓
- [x] `DiffusionExperiment` (Stejskal–Tanner, D + rH) so `diffusion1d` has a replacement
- [x] legacy `relaxation1d.jl` / `tract.jl` / `calibration1d.jl` / `diffusion1d.jl` /
      `regions1d.jl` removed; the new implementations take those names
- [x] registered with the analysis-dispatch registry (R1/R2 relaxation, nutation
      calibration), replacing the registrations the legacy files made
- [x] `pickregion` — region selection on its own, overlaying every supplied spectrum —
      and `Exchange1D._prompt_integration!` now uses it instead of readline prompts,
      feeding the existing `integrate!(prob, peakppm, noiseppm, ppmwidth)` unchanged
- [ ] dispatch rules for TRACT (needs a TROSY/anti pair), diffusion (needs the gradient
      list) and STD/kinetics (their default single region is registrable in principle -
      revisit)

Phase 3.5 — **first real-session feedback round**  ← THIS ITERATION
- [x] noise is now a single draggable position (vline + narrow drag target), not an
      independently-sized region: the noise region always matches the width of whichever
      signal region is being reduced (matching widths is what makes it a direct estimate
      of that region's integration noise) — `Dataset1D.noisecenter::Float64` replaces the
      old `noise::Region` field
- [x] error bars corrected to standard practice: std, **across planes**, of the noise
      region's integral (was: per-point noise inside one trace, scaled by √N) — see
      `reduce_region`
- [x] intensities scaled down by that noise level (every point then has σ = 1 by
      construction), so numbers read in signal-to-noise units instead of raw integrated
      intensities
- [x] default region width fixed at 0.05 ppm, centred on the tallest peak
      (`defaultregion`) — was defaulting to the full spectral span
- [x] outer-layout `colsize!`/`rowsize!` (`Auto(false, ratio)`, matching GUI2D) so the
      window actually fills on resize
- [x] fitting toggle is now a `Toggle`, not a `Button`
- [x] region dropdown only shown when there is more than one region
- [x] interactive region add/rename/delete, matching the 2D fitting GUI's key bindings:
      `A` (tap or press-drag-release) adds, `R` renames, `D` deletes, plus matching
      buttons and on-screen help text; `std1d`/`kinetics1d` no longer require `regions`
      up front
- [x] Left/Right arrow keys step through spectra (renamed from "planes" throughout the UI)
- [x] active-region highlighting on hover (and on click/drag), synced with the dropdown
- [x] moving/resizing/switching the active region refreshes the fit-panel axis limits
- [x] multi-series result plotting (`ResultSeries`): TRACT's TROSY/anti-TROSY and STD's
      saturation frequencies each get their own colour, matched between points and fit
      line, with an inline coloured label (no formal `Legend` yet)
- [x] tidied `summary_text` (aligned, rounded, units) — still plain text; **open
      question, not resolved here**: what should the canonical output format be across
      1D *and* 2D (text vs CSV, units, precision)? Needs a decision, then one shared
      formatter.
- [x] fixed the `Exchange1D` crash: `traces_from_spec` assumed a 2D (shift × plane)
      array and indexed `Y[:, i]`, which failed on genuinely N-D pseudo-3D/-nD Exchange1D
      spectra (e.g. off-resonance R1ρ, shift × spinlock × delay); it now reshapes any
      array whose first dimension is the chemical shift into a flat list of planes

Phase 4 — Tier-2 exchange second window (`ExchangeModel`, parameter table with fix/global)
as a front-end onto `ExchangeProblem`/`FitResult`; the reduced dispersion / Z-profiles in
the Tier-1 result panel.

Phase 5 — additional reductions (NMF kinetics), qNMR/PULCON, temperature calibration;
lift a shared `AnalysisCore` out of GUI2D + Analysis1D; a real `Legend` for multi-series
result plots; per-region-edge drag-resize (currently: drag recentres, textbox resizes).

Phase 3.6 — **restructure onto the GUI2D pattern**  ✓
- [x] one self-contained `expt-<name>.jl` per experiment (entry point, type, interface,
      science, presentation), included from `experiments.jl`; `loaders.jl` reduced to the
      NMRData adapters in `nmrdata.jl`, and each experiment's acquisition-parameter
      physics moved beside the science that uses it
- [x] `RegionResult` replaces `SeriesResult`: the uniform `parameters`/`postparameters`
      pair GUI2D's `Peak` carries, filled by `postfit!`/`postfitglobal!`, with
      `primaryparam` naming the headline quantity. `analyse` returns a
      `Vector{RegionResult}` for every experiment, which also removes the Observable
      element-type fragility the old `(; series, summary)` shape caused
- [x] the summary formatter is now generic over those two dictionaries — the six
      per-experiment `_summary_extra` methods are gone — with one `PARAM_UNITS` /
      `PARAM_LABELS` table. **This is the first step on the open question below**: the
      remaining work is to merge it with `gui2d/summary.jl`'s `PARAM_LABELS`
- [x] STD no longer overrides `analyse`; it contrasts and fits in `postfitglobal!`, the
      way `cest2d` does in `postfit!`, and its four presentation overrides are deleted
- [x] `ResultVisualisation` trait (`completeresultstate!` / `resultpanel!` /
      `plotresult!`), mirroring GUI2D's `VisualisationStrategy`, so a future experiment
      needing a different panel — CEST Z-profiles, the exchange window — has somewhere to
      go. `gui.jl` still contains no per-experiment special-casing at all
- [x] regions have one representation (`Region`) rather than two
- [x] `files.jl`: `results.csv` plus a region round-trip, as in 2D, so a multi-region STD
      or kinetics session is reproducible; `experimentinfo` supplies its header
- [x] internals renamed to house style (no underscores, no `_` prefixes)
- [x] `docs/src/advanced/creating_1d_analyses.md`

Still open after this iteration:
- **Canonical output format across 1D and 2D** (text vs CSV, units, precision). 1D now has
  a `results.csv` following 2D's conventions and a single units/labels table, but the two
  tables have not been merged and the question itself is not settled.
- The region dropdown is created only when the experiment starts with more than one
  region, so regions added later (with `A`, or by loading a saved list) are selectable by
  clicking the plot but do not appear in a menu.
