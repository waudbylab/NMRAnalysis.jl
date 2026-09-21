# Automatic Analysis

`analyse()` looks at the annotations embedded in one or more experiment files and either
runs the matching analysis directly, or offers a menu when more than one could apply. It
saves you from having to know in advance which function handles a given experiment —
useful when working through a folder of mixed experiment types, or when scripting a
pipeline that shouldn't need to know what's coming next.

!!! note "1D only, for now"
    Automatic dispatch currently covers 1D experiments only. 2D analyses (`fit2d()`,
    `relaxation2d()`, and the rest) are always called directly — see [2D
    Analysis](2d/overview.md). Dispatch for 2D experiments is planned.

```julia
using NMRAnalysis

# Single file
result = analyse("data/101")

# Multiple files — e.g. a CEST series plus an R1 experiment for a joint exchange fit
results = analyse(["data/101", "data/102", "data/103"])
```

## How it works

1. **Classification**: each input file is classified by its `types` (e.g. `"1d"`,
   `"calibration"`, `"r1rho"`) and `features` (e.g. `"R1"`, `"nutation"`,
   `"on_resonance"`), read from annotations embedded in the pulse sequence.
2. **Matching**: registered analysis rules are checked against the classified files.
3. **Selection**: if exactly one rule matches, it runs immediately; if several match, an
   interactive menu lets you choose.
4. **Execution**: the selected analysis runs and its result is returned.

Because matching depends on annotations carried by the pulse sequence, `analyse()` only
recognises experiment types whose sequences include them — it isn't inferring the
experiment type from the raw data itself.

Dispatch is a convenience, not a requirement: every routine can be called directly, and
none of them needs annotations to run. Where an annotation is missing, the routine asks
for what it needs instead. See
[Where the parameters come from](1d/overview.md#Where-the-parameters-come-from).

## What it currently detects

| Input | Runs |
|---|---|
| One or more 1D nutation calibration experiments | [Pulse Calibration](1d/calibration.md) |
| A 1D R1 or R2 relaxation experiment | [Relaxation](1d/relaxation.md) |
| One or more on-resonance R1ρ experiments | [R1ρ Relaxation Dispersion](r1rho.md) |
| CEST and/or off-resonance R1ρ experiments, together with any accompanying R1 or on-resonance R1ρ files, and any nutation calibrations | [Exchange Fitting](exchange1d/overview.md) |

Several nutation calibrations selected together are analysed as one calibration curve
rather than one at a time. Where they are selected alongside CEST or R1ρ experiments they
are offered both ways: on their own as a calibration to inspect, and as part of the
exchange analysis, where they are not fitted as data but used for the spin-lock field
strengths and the B₁ inhomogeneity. Taken that second way they are fitted without a window,
and their own results are written into a `calibration/` folder inside the exchange fit's
output when that fit is saved, so there is still something to check them by. Without them
the exchange fit falls back to each experiment's own reference pulse and an assumed 5%
inhomogeneity.

Diffusion and TRACT don't currently register with `analyse()` — call
[`diffusion1d()`](1d/diffusion.md) or [`tract()`](1d/tract.md) directly.

## Extending it

`analyse()` is extensible: analysis modules register themselves at load time, so any
package built on this mechanism is automatically discovered without changes here. See
[Analysis Rules](../advanced/analysis_rules.md) for how to register a new analysis type.
