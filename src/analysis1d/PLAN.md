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

- **Global post-fit** (`postprocess`) — combine series into a derived result
  (TRACT τc from ΔR₂; STD epitope normalisation). Returns a scalar/table shown inline,
  or — where it is an interactive model fit (exchange) — drives a second window.
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
out of both (start parallel, refactor toward shared core once 1D is proven).

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
- [x] grouping + curve-fit pipeline (noise-weighted) → `SeriesResult`
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

Phase 3 — **replace the readline routines + hook into dispatch**  ← THIS ITERATION
- [x] `DiffusionExperiment` (Stejskal–Tanner, D + rH) so `diffusion1d` has a replacement
- [x] legacy `relaxation1d.jl` / `tract.jl` / `calibration1d.jl` / `diffusion1d.jl` /
      `regions1d.jl` removed; the new implementations take those names
- [x] registered with the analysis-dispatch registry (R1/R2 relaxation, nutation
      calibration), replacing the registrations the legacy files made
- [x] `pickregion` — region selection on its own, overlaying every supplied spectrum —
      and `Exchange1D._prompt_integration!` now uses it instead of readline prompts,
      feeding the existing `integrate!(prob, peakppm, noiseppm, ppmwidth)` unchanged
- [ ] dispatch rules for TRACT (needs a TROSY/anti pair), diffusion (needs the gradient
      list) and STD/kinetics (need named regions)

Phase 4 — Tier-2 exchange second window (`ExchangeModel`, parameter table with fix/global)
as a front-end onto `ExchangeProblem`/`FitResult`; the reduced dispersion / Z-profiles in
the Tier-1 result panel.

Phase 5 — additional reductions (NMF kinetics), qNMR/PULCON, temperature calibration;
lift a shared `AnalysisCore` out of GUI2D + Analysis1D.
