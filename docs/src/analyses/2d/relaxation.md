# Relaxation Analysis (T1 / T2)

The `relaxation2d` function measures R1 or R2 relaxation rates from a series of 2D
spectra recorded with increasing relaxation delays. Peak amplitudes are fitted to a
mono-exponential decay:

```math
I(\tau) = A \exp\!\left(-R\tau\right)
```

where ``R`` is the relaxation rate (s⁻¹) and ``A`` is the peak amplitude. The software
does not distinguish between R1 and R2 — the appropriate interpretation depends on the
experiment used to collect the data.

![Screenshot from relaxation analysis](../../assets/screenshot-relaxation2d.png)

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

Clicking **Save to folder** writes all results to `results.csv`. Alongside peak
positions, linewidths and the amplitude for each delay, the derived columns are:

| Column | Description |
|--------|-------------|
| `R`, `R_err` | Fitted relaxation rate R (s⁻¹) and uncertainty |
| `A`, `A_err` | Fitted amplitude A and uncertainty |

The rate is labelled generically as `R`; the software does not distinguish R₁
from R₂. See [Peak Lists and Output Files](peaklistformats.md) for the full format.

Plot R against residue number with [`summaryplot`](summary.md). Pass an appropriate
`ylabel` to label the axis for your specific experiment:

```julia
# T2 / R2 measurement
fig = summaryplot(expt; ylabel="R₂ / s⁻¹")
fig = summaryplot("results/"; param=:R, ylabel="R₂ / s⁻¹")

# T1 / R1 measurement
fig = summaryplot(expt; ylabel="R₁ / s⁻¹")
```

## Noise Estimation

Peak amplitude uncertainties are estimated from the scatter of the spectral noise across
the series of experiments.
