# Simple fitting

`fit2d` fits peaks in one 2D spectrum, or a series, with no model: it reports each peak's
position, linewidths and amplitudes. Use it to measure shifts and linewidths, to build a
peak list for another analysis, or to compare intensities across a series.

Each peak is a 2D Lorentzian, with the window function used in processing applied.
Amplitudes are divided by the number of scans, so they compare across spectra recorded
with different numbers of scans.

## Usage

```julia
using NMRAnalysis
fit2d(11)          # one spectrum
fit2d(11:13)       # a series
```

### Example

The ¹H,¹⁵N HSQC of ubiquitin:

```julia
fit2d(example("2d_hn_ubq"))
```

```@setup fit2d
using NMRAnalysis
NMRAnalysis.screenshot("fit2d.png") do
    fit2d(example("2d_hn_ubq"))
end
```

![The fit2d window, before any peaks are picked](fit2d.png)

## Output

**Save** writes these parameters for each peak to `results.csv`, a plot of each peak's fit
to `peaks/`, and a summary plot against residue number to `summary.pdf`. See
[Output](overview.md#Output) for the other files.

| Parameter | Meaning |
|---|---|
| `x`, `y` | Position, ppm |
| `R2x`, `R2y` | Linewidths, s⁻¹ |
| `amp` | Amplitude in each spectrum |

The amplitudes, one per spectrum, are in `series.csv`. To fit a model to them, use
[`relaxation2d`](relaxation.md), [`recovery2d`](magnetisationrecovery.md) or
[`modelfit2d`](modelfit.md).
