# Magnetisation recovery

`recovery2d` fits the amplitude of each peak in an inversion or saturation recovery
experiment to

```math
I(T) = A \left[ 1 - C \exp(-RT) \right],
```

where ``R`` is the recovery rate and ``C`` is 2 for a complete inversion and 1 for
saturation.

![The recovery2d window](../../assets/recovery2d.png)

## Usage

```julia
using NMRAnalysis
recovery2d("33"; T=[0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
recovery2d("33"; T="vdlist.txt")
```

The delays are found as for [`relaxation2d`](relaxation.md). The window opens on the plane
with the longest delay, where the peaks are easiest to pick.

## Output

**Save** writes these parameters for each peak to `results.csv`, a plot of each peak's fit
to `peaks/`, and a summary plot against residue number to `summary.pdf`. See
[Output](overview.md#Output) for the other files.

| Parameter | Meaning |
|---|---|
| `R` | Recovery rate, s⁻¹ |
| `A` | Equilibrium amplitude |
| `C` | Recovery factor |
