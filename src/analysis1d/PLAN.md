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
2. **Series model** — quantity vs evolution parameter → derived parameters: v1 is
   curve-fit only (exponential, recovery, damped sinusoid, …), continuous fit-axis. A
   *contrast* shape (categorical slices, e.g. a reference-vs-saturated fraction — the 1D
   analogue of the existing `hetnoe2d` reference/saturated experiment) was prototyped for
   STD and removed with it (see "Removed" below); revisit if a future experiment needs it.
3. **Visualisation strategy** — usually derived from the series model; overridable.

Plus:

- **Post-fit** (`postfit!` / `postfitglobal!`) — derive further quantities from a fit and
  record them on the result (TRACT τc from ΔR₂), the same hooks GUI2D uses.
  `postfitglobal!` is the one to use when the derivation spans several series.
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
| kinetics | N named | `time, run` | `time` | `run` | Integrate | NoFitting (v1) | – |

Entry points take the names of the routines they replace (`relaxation1d`, `tract`,
`calibration1d`, `diffusion1d`, plus `kinetics1d`). Each opens the GUI; passing
`integration = (; peakppm, noiseppm, ppmwidth)` skips it and analyses that region
directly, so a chosen region can be replayed from a script. The triple is deliberately
the same one `Exchange1D` stores as `prob.integration`.

## Removed: STD

An STD experiment (`std1d`) was built and shipped in earlier phases of this iteration,
then removed: its entry point, `std1d(spec, sat, tsat; …)`, assumed the entire experiment
— every saturation frequency crossed with every saturation time — already existed as
planes within *one* loaded file, with the caller supplying flat `sat`/`tsat` vectors
matching that file's flattened plane order. Real STD acquisition is normally the other
shape: a separate pseudo-2D dataset **per saturation time**, each one internally arraying
saturation frequency (usually on/off-resonance, sometimes a handful more) via an fq-list —
structurally the same as TRACT's `tract(trosy, antitrosy)`, which takes two separate files
and combines them itself, not one experiment's worth of `std1d`. Getting this right needs
proper design (how many files, what's arrayed within each, whether saturation time is ever
itself arrayed in one file) before reintroducing it, rather than patching the existing
entry point's assumptions. The analysis-side pipeline it used — reduce → group by
saturation frequency → contrast against a reference → fit a buildup curve → normalise into
an epitope map, via `postfitglobal!` (see Phase 3.6/3.7 below) — is still a reasonable
target shape for whatever `std1d` becomes; it is the *loading* step that needs rethinking.

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
      (later removed entirely - see Phase 3.6; hovering selects a region, so the dropdown
      was redundant)
- [x] interactive region add/rename/delete, matching the 2D fitting GUI's key bindings:
      `A` (tap or press-drag-release) adds, `R` renames, `D` deletes, plus matching
      buttons and on-screen help text; `std1d`/`kinetics1d` no longer require `regions`
      up front
- [x] Left/Right arrow keys step through spectra (renamed from "planes" throughout the UI)
- [x] active-region highlighting on hover (and on click/drag)
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
      per-experiment `_summary_extra` methods are gone. (Superseded by Phase 3.9 below:
      `PARAM_UNITS`/`PARAM_LABELS` briefly grew into one flat table for every experiment's
      parameters, which was wrong for the same reason `_summary_extra` was — corrected to
      one small shared table plus a per-experiment override table each.)
- [x] STD no longer overrides `analyse`; it contrasts and fits in `postfitglobal!`, the
      way `cest2d` does in `postfit!`, and its four presentation overrides are deleted
- [x] `ResultVisualisation` trait (`completeresultstate!` / `resultpanel!` /
      `plotresult!`), mirroring GUI2D's `VisualisationStrategy`, so a future experiment
      needing a different panel — CEST Z-profiles, the exchange window — has somewhere to
      go. `gui.jl` still contains no per-experiment special-casing at all
- [x] regions have one representation (`Region`) rather than two
- [x] `files.jl`: `results.csv` plus a region round-trip, as in 2D, so a multi-region STD
      or kinetics session is reproducible; `experimentinfo` supplies its header
- [x] region dropdown removed: hovering a region makes it active, so the menu was
      redundant - and it was only ever created when the experiment *started* with more
      than one region, so it never listed regions added later
- [x] internals renamed to house style (no underscores, no `_` prefixes)
- [x] `docs/src/advanced/creating_1d_analyses.md`

Phase 3.7 — **second feedback round**  ✓
- [x] TRACT's `postfitglobal!` recorded η/τc on both series of the TROSY/anti-TROSY pair,
      so the results panel - which prints every series of the active region - showed the
      derived summary twice; recorded once, on the TROSY member, instead
- [x] the redundant region dropdown removed: hovering (or clicking/dragging) a region
      already makes it active, which is how the dropdown's own selection was driven in the
      first place, and it was only ever populated from the region list at GUI-open time,
      so a region added afterwards was never listed there anyway
- [x] `pickregion`'s noise band now tracks the integration width the way `gui!`'s does
      (matching widths is what makes the noise integral a direct estimate of the signal
      region's noise), floored at the default region width, rather than a fixed size; and
      peak position, noise position and width each get a two-way-bound text box
      (`ppmbox!`) for exact entry alongside dragging
- [x] the results panel's group headers (`resultsheader`) and derived-quantity headers
      (`secondarytext`, via the new `derivedheader` hook) are now bold `RichText`
      (`boldheader`), so the panel reads as labelled sections rather than one
      undifferentiated block of text - `state[:resultspanel]` (see below) carries
      `RichText` throughout, including its blank states, for the same
      Observable-element-type reason `RegionResult` exists
- [x] regions added with `A` default to `peak1`, `peak2`, … rather than `region1`,
      `region2`, … - so the panel's "Region: peak2" heading doesn't repeat "region"
      (`defaultregion`'s own "signal" default, used by every experiment's *first* region,
      is untouched)
- [x] TRACT's default region is now the conventional 7.5-9.5 ppm amide envelope
      (`defaultamideregion`), not a narrow window on the tallest peak - TRACT integrates
      the bulk backbone signal, and a single-peak default made little scientific sense
      for it specifically
- [x] shift+scroll resizes the region under the cursor (about its own centre), in both
      `gui!`'s main window and `pickregion`'s standalone popup - the same operation the
      width textbox performs, via the wheel; plain scroll is untouched, so the axis's own
      x-zoom still works normally
- [x] the spectrum axis title reads "Spectrum: 0.4 s delay" rather than a bare
      "0.4 s delay"

Phase 3.8 — **third feedback round: one info panel, real alignment**  ✓
- [x] the right-hand info panel (which region, its bounds, the fit, and the quantities
      derived from it) is now `state[:resultspanel]`, a single `RichText` Label rather
      than three separately-sized ones. Each Label in a GridLayout column auto-sizes its
      row from its own reported height, and that reporting does not reliably track a
      `RichText`'s actual rendered height as its content changes - so a later row could
      start before an earlier one had actually finished, overlapping it (TRACT's η/τc
      block visibly colliding with its own trosy/anti blocks). One Label sidesteps the
      question entirely: everything after the heading is just more content appended to
      the same growing block, with nothing else positioned relative to it to collide with
- [x] `panelwidth` computes the label-column width once, across every block shown for the
      active region, rather than each block aligning only to its own labels - so
      "Amplitude" in a TROSY block and "Correlation time (τc)" in TRACT's derived block
      now share one column
- [x] `paramlabel`/`paramunit` are experiment-dispatched, not one global table: `:A` is
      safely "Amplitude" everywhere, but `:R` is deliberately *not* given a global
      "Relaxation rate" label, since nutation's damped-sinusoid model uses the same bare
      symbol for a decay rate that is not conventionally called that - relaxation and
      TRACT each override it themselves, where it is actually true
- [x] fixed a real bug this surfaced: TRACT's derived (η/τc) block was headed "trosy",
      because `postfitglobal!` has to pick one series of the pair to record them on and
      the header defaulted to that series' own group name. New `derivedheader` hook,
      overridden per experiment, decouples "whose `RegionResult` holds this" from
      "what to call it" - TRACT now says "TRACT results" regardless of which member
      carries the values
- [x] `groupheader` (new hook, default = `groupname`) lets TRACT show "TROSY"/"Anti-TROSY"
      in the panel, matching the wording `seriesnames` already uses for the plot legend,
      instead of the raw `:trosy`/`:anti` group values
- [x] every panel heading is now colon-suffixed and the "=" between name and value is
      gone, matching plain name/value columns rather than an equation
- [x] `word_wrap` removed from the results-panel Label - it made the Label's wrap width
      follow the cell GridLayout had given it while that cell's own size came from the
      Label's *reported* size, a feedback loop that does not reliably settle for
      reactively-changing content; symptoms were a large fixed gap that didn't track the
      window height and text reading as centred rather than left-aligned
- [x] the same row's *height*, however, turned out not to be reliably reported by Makie's
      own RichText boundingbox/autosize machinery either - fixed for a short,
      single-block panel (`relaxation1d`) but still stuck too short for a longer one with
      a derived-quantity block too (`calibration1d`). Rather than continue chasing that
      machinery, the row height is now computed directly (`PANEL_LINE_HEIGHT` ×
      newline-count in the flattened text) and set explicitly via `rowsize!`, so it no
      longer depends on Makie's auto-sizing for this element at all
- [x] shift+scroll resizes `state[:active][]` directly rather than requiring the cursor
      to stay precisely over the region's own span - more robust for a physical mouse
      wheel's tendency to nudge the cursor slightly as it turns
- [x] four more fit parameters given display names (`C` "Recovery factor", `D` "Diffusion
      coefficient", `ν` "Nutation frequency", `k` "Buildup rate") - added to the shared
      `PARAM_LABELS` table at the time; corrected in Phase 3.9 below, since none of the
      four is actually shared

Phase 3.9 — **fourth feedback round: parameter names belong with their experiment**  ✓
- [x] every experiment-specific entry moved out of the shared `PARAM_UNITS`/
      `PARAM_LABELS` tables in `visualisation.jl` and into a `<EXPT>_PARAM_LABELS`/
      `<EXPT>_PARAM_UNITS` pair in that experiment's own file, with `paramlabel`/
      `paramunit` overridden there to check it first. Prompted by nutation's "B₁ inhom."
      label turning up nowhere in `expt-nutation.jl` - the override mechanism built for
      `:R` (Phase 3.7) was right, but every symbol added after it (TRACT's η/τc,
      nutation's pulse90/inhomogeneity, diffusion's D/rH/viscosity, STD's whole
      buildup/epitope vocabulary) went into the shared table instead of using it, the
      same "everything about an experiment lives with that experiment" violation the
      whole `expt-*.jl` split exists to prevent. The shared tables now hold only what is
      genuinely cross-experiment: `:A`'s label (the same "Amplitude" everywhere it's
      fitted) and `:R`'s unit (a rate is a rate regardless of which experiment fits one -
      only its *label* differs by experiment, which is why that stays overridden per
      experiment rather than joining `:A` here)
- [x] `RecoveryModel` moved from `seriesmodels.jl` (shared-model file) to
      `expt-relaxation.jl`: it was never actually shared, only `RelaxationExperiment`'s
      `ir=true` path ever constructs one - the same misplacement as the labels, just for
      a model instead of a name

Still open after this iteration:
- **Canonical output format across 1D and 2D** (text vs CSV, units, precision). 1D now has
  a `results.csv` following 2D's conventions and a small shared units/labels table (plus
  one per experiment for what isn't shared), but the 1D and 2D tables have not been
  merged and the question itself is not settled.
