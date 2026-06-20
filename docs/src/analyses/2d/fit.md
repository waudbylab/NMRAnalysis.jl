# 2D Peak Fitting

```@docs
fit2d
```

The `fit2d` function opens an interactive GUI for fitting peaks in a single 2D spectrum
or a series of 2D spectra. Each peak is fitted to a 2D Lorentzian lineshape; no physical
model is applied to the amplitudes across spectra — all measured quantities (peak
positions, linewidths, and amplitudes) are reported directly.

This makes `fit2d` a flexible starting point for many different experiments. Typical use
cases include:

- **Single spectrum**: measure chemical shifts and linewidths, or build a peak list for
  use with other analysis functions
- **Intensity series**: compare peak amplitudes across a concentration titration,
  temperature series, or time course without assuming a model
- **Saturation transfer difference**: compare reference and irradiated spectra
- **Solvent exposure**: measure intensity ratios between samples with and without a
  paramagnetic probe (see also [`pre2d`](@ref) for a full PRE analysis)

For data that follow a known physical model, consider the more specific functions:
[`relaxation2d`](@ref) (exponential decay), [`recovery2d`](@ref) (magnetisation
recovery), or [`modelfit2d`](@ref) (arbitrary equation).

![Screenshot from peak fitting to a single 2D](../../assets/intensity2d.png)

## Usage

```julia
using NMRAnalysis

# Single spectrum
fit2d("expno/pdata/1")

# Series of spectra — track peaks across multiple experiments
fit2d([
    "expno1/pdata/1",
    "expno2/pdata/1",
    "expno3/pdata/1",
])
```

The input path(s) should point to processed Bruker data directories (i.e. `pdata/N`
subdirectories containing `2rr` files).

## Output

Clicking **Save to folder** writes peak positions, linewidths, and amplitudes to
`fit.peaks`. Each row corresponds to one peak; columns include:

| Column | Description |
|--------|-------------|
| `label` | Peak label (editable in the GUI) |
| `residue` | Residue number extracted from label |
| `x`, `y` | Fitted chemical shifts (ppm) |
| `R2x`, `R2y` | Fitted linewidths in the direct and indirect dimensions (s⁻¹) |
| `amp1`, `amp2`, … | Fitted amplitude for each input spectrum |

For a single-spectrum analysis, there is one amplitude column. For a series, there is one
column per spectrum.
