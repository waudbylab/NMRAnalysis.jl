# Automatic Analysis

`analyse()` looks at the annotations embedded in one or more experiment files and either
runs the matching analysis directly, or offers a menu when more than one could apply. It
saves you from having to know in advance which function handles a given experiment —
useful when working through a folder of mixed experiment types, or when scripting a
pipeline that shouldn't need to know what's coming next.

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

## What it currently detects

| Input | Runs |
|---|---|
| A 1D nutation calibration experiment | [Calibration](calibration.md) |
| A 1D R1 or R2 relaxation experiment | [Relaxation](relaxation.md) |
| One or more on-resonance R1ρ experiments | [R1ρ Relaxation Dispersion](r1rho.md) |
| CEST and/or off-resonance R1ρ experiments, together with any accompanying R1 or on-resonance R1ρ files | [Exchange Fitting](exchange1d.md) |

Diffusion and TRACT don't currently register with `analyse()` — call
[`diffusion1d()`](diffusion.md) or [`tract()`](tract.md) directly.

## Extending it

`analyse()` is extensible: analysis modules register themselves at load time, so any
package built on this mechanism is automatically discovered without changes here. See
[Analysis Rules](../advanced/analysis_rules.md) for how to register a new analysis type.
