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

![Screenshot of relaxation fitting](../../assets/relaxation-demo.mov)

## Usage

```julia
using NMRAnalysis

# Pseudo3D data (experiment 20) whose delays are in its vdlist or annotations
relaxation2d("20")

# Pseudo3D data with a list of relaxation times (in seconds)
relaxation2d("20"; relaxationtimes="20/relaxation-times.txt")

# Or provide a list of 2D planes and associated relaxation times
relaxation2d(["11", "12", "13", "14", "15"];
             relaxationtimes=[0.010, 0.030, 0.060, 0.100, 0.200])

# A delay counted in loops of a vclist, each loop 16 ms long
relaxation2d("21"; cycletime=0.016)

# Or loop counts given directly
relaxation2d("21"; ncyc=[0, 2, 4, 8, 16], cycletime=0.016)
```

The number of input planes must match the number of relaxation delays. The delays are
looked for in this order:

1. `relaxationtimes`, a vector or the path of a file holding one per line
2. `ncyc` multiplied by `cycletime`, for delays counted in loops
3. the pulse-sequence annotation `relaxation.duration`
4. the `vdlist`
5. the `vclist` multiplied by `cycletime`
6. a question, before the window opens

Where a loop count is used and `cycletime` isn't given, you are asked for it. The second
argument may still be given positionally, as `relaxation2d("20", delays)`.

## Excluding planes from the fit

If one or more planes in the series should not contribute to the fitted rate — pass their 1-based indices via `skipplanes`:

```julia
relaxation2d(files; relaxationtimes=delays, skipplanes=[1, 5])
```

All spectra are still loaded and displayed. Skipped planes appear as open grey markers
in the peak-fit plot and are labelled **[skipped]** in the slice title; they are not
used when fitting R or A. The full delay list, including times for the skipped planes,
must always be supplied.

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
fig = summaryplot("results/"; param=:R, ylabel="R₂ / s⁻¹")

# T1 / R1 measurement
fig = summaryplot("results/"; param=:R, ylabel="R₁ / s⁻¹")
```

## Noise Estimation

Peak amplitude uncertainties are estimated from the scatter of the spectral noise across
the series of experiments.
