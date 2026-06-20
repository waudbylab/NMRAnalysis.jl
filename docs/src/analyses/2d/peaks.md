# Peak Analysis

```@docs
peaks2d
```

The `peaks2d` function opens an interactive GUI for picking peaks in a single 2D spectrum
or a series of 2D spectra, fitting each peak to a 2D Lorentzian lineshape. No physical
model is applied to the amplitudes across spectra — all measured quantities (peak
positions, linewidths, and amplitudes) are reported directly.

Typical use cases include:

- **Chemical shift mapping**: track peak positions across a titration or temperature series
- **Linewidth analysis**: measure R2x and R2y linewidths across conditions
- **Intensity ratios**: compare peak amplitudes without assuming a model (e.g. saturation
  transfer difference, solvent exposure)
- **Reference peak list**: build a peak list for use with other analysis functions

![Screenshot from peak fitting to a single 2D](../../assets/intensity2d.png)

## Usage

```julia
using NMRAnalysis

# Single spectrum
peaks2d("expno/pdata/1")

# Series of spectra — track peaks across multiple experiments
peaks2d([
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
