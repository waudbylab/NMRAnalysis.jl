# Relaxation Analysis (T1 / T2)

```@docs
relaxation2d
```

The `relaxation2d` function measures R1 or R2 relaxation rates from a series of 2D
spectra recorded with increasing relaxation delays. Peak amplitudes are fitted to a
mono-exponential decay:

```math
I(\tau) = A \exp\!\left(-R\tau\right)
```

where ``R`` is the relaxation rate (s⁻¹) and ``A`` is the peak amplitude. The software
does not distinguish between R1 and R2 — the appropriate interpretation depends on the
experiment used to collect the data.

<!-- screenshot: assets/relaxation2d.png -->

## Usage

```julia
using NMRAnalysis

# Inline relaxation delays (in seconds)
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1", "15/pdata/1"],
    [0.010, 0.030, 0.060, 0.100, 0.200]
)

# Read delays from a text file (one value per line; lines beginning with # are ignored)
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1"],
    "vclist.txt"
)
```

The number of input spectra must match the number of relaxation delays.

## Output

Clicking **Save to folder** writes fitted relaxation rates to `fit-results.txt`. Columns
include:

| Column | Description |
|--------|-------------|
| `label` | Peak label |
| `R_value` | Fitted relaxation rate R (s⁻¹) |
| `R_uncertainty` | Uncertainty in R (s⁻¹) |
| `A_value` | Fitted amplitude A |
| `A_uncertainty` | Uncertainty in A |

The `fit.peaks` file additionally contains fitted peak positions, linewidths, and the
amplitude for each delay.

## Noise Estimation

Peak amplitude uncertainties are estimated from the scatter of the spectral noise across
the series of experiments.
