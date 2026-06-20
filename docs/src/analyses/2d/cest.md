# CEST

```@docs
cest2d
```

CEST (Chemical Exchange Saturation Transfer) experiments detect the presence of minor
conformational states through their effect on the major-state peak intensity. A weak
radiofrequency field is applied at a variable saturation offset; when the offset matches
the resonance frequency of a minor-state spin, saturation is transferred to the
major-state peak via chemical exchange, reducing its intensity. The resulting
Z-spectrum (also called a CEST profile) reveals both the major-state and minor-state
chemical shifts.

<!-- screenshot: assets/cest2d.png -->

## Usage

The input is a single pseudo-3D dataset where the first plane is the reference spectrum
(recorded without saturation) and subsequent planes are saturation spectra recorded at
increasing saturation offsets. The saturation frequencies are read automatically from the
`fq3list` file in the experiment directory.

```julia
using NMRAnalysis

cest2d("11/pdata/1"; B1=15, Tsat=0.3)
```

- `B1`: Saturation field strength in Hz. Typical values are 5–50 Hz for ¹⁵N CEST.
- `Tsat`: Saturation time in seconds.

## Output

Clicking **Save to folder** writes all results to `results.csv`. Each row is one
peak, with the per-offset amplitudes (`amp[1]`, `amp[2]`, …) that make up the
Z-spectrum (normalised intensity ``I(\omega_\text{sat})/I_0``), plus the fitted
`R1`/`R2` rates and uncertainties. See
[Peak Lists and Output Files](fileformats.md) for the full format.
