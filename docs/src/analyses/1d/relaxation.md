# Relaxation (R₁, R₂)

`relaxation1d` fits a series of 1D spectra recorded with increasing relaxation delay, and
returns the relaxation rate. It handles both a plain decay, as in a T₂ or a saturation
experiment, and an inversion recovery.

## Running the analysis

```julia
using NMRAnalysis

results = relaxation1d("11")
```

That reads experiment 11, works out the delays and the model (see below), and opens the
[analysis window](overview.md#Using-the-analysis-window). Put the region over your signal,
put the noise marker somewhere empty, check the fit, and press **Save**.

<!-- TODO screenshot (2): docs/src/assets/relaxation1d-window.png
     relaxation1d on a ¹⁹F R₁ inversion recovery: the region over the peak, the recovery
     curve fitted below, and the results panel showing the rate. Replaces the old
     R2-fit.png / IR-fit.png, which show the Plots.jl output of the retired routine. -->

## Delays and model

The delays and the model are found in the usual
[order of precedence](overview.md#Where-the-parameters-come-from). Nothing has to be
annotated: if the delays cannot be found you are asked for them, and if the model cannot
be found you are offered the choice.

| | From an argument | From an annotation | From the acquisition parameters | Otherwise |
|---|---|---|---|---|
| Delays | `tau=[0.01, 0.05, …]` | `relaxation.duration` | the `vdlist` | you are asked |
| Model | `model=:exponential` or `:recovery` | `relaxation.model` | | you choose from a menu |

```julia
relaxation1d("11"; tau=[0.01, 0.03, 0.06, 0.1])   # delays given explicitly
relaxation1d("11"; model=:recovery)               # fit a recovery curve
relaxation1d("11"; tau=collect(0.01:0.01:0.1), model=:exponential)   # nothing asked
```

Delays are always in seconds.

## What is fitted

An exponential decay,

```math
I(t) = A \exp(-Rt),
```

or, for `model = :recovery`, an inversion or saturation recovery,

```math
I(t) = A \left[ 1 - C \exp(-Rt) \right].
```

| Parameter | Meaning |
|---|---|
| `A` | Amplitude |
| `R` | Relaxation rate, s⁻¹ |
| `C` | Recovery factor, for a recovery fit only. Ideally 2 for a complete inversion |

The relaxation time is `1/R`. `param(results[1], :R)` gives the rate with its uncertainty,
so `1 / param(results[1], :R)` gives the time with the uncertainty propagated.

## Automatic dispatch

A pulse sequence annotated with types `"1d"` and `"relaxation"`, and a feature of `"R1"` or
`"R2"`, is recognised by [`analyse`](../analyse.md):

```julia
analyse("11")
```

## Notes

For an inversion recovery, step to the last spectrum before positioning the region: that
is where the signal has recovered most and is easiest to see.

The uncertainty on each integral comes from the noise marker, so it depends on a flat
baseline there. See [Uncertainties](overview.md#Uncertainties).
