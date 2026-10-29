# CEST

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
increasing saturation offsets.

```julia
using NMRAnalysis

cest2d("11"; B1=15, Tsat=0.3)

# Read everything from the pulse-sequence annotations
cest2d("11")
```

- `B1`: Saturation field strength in Hz. Typical values are 5–50 Hz for ¹⁵N CEST. If you
  don't give it, it is calculated from the annotated saturation power (`cest.power`),
  calibrated against the reference pulse of the annotated channel (`cest.channel`).
- `Tsat`: Saturation time in seconds, otherwise read from the `cest.duration` annotation.
- `offsets`: Saturation offsets in ppm, one per plane, otherwise read from the
  `cest.offset` annotation or the `fq3list`.
- `skipplanes`: Saturation planes to leave out of the fit. The reference, plane 1, can't
  be skipped.

Anything that can't be found is asked for before the window opens.

## Output

Clicking **Save to folder** writes all results to `results.csv`. Each row is one
peak, with the per-offset amplitudes (`amp[1]`, `amp[2]`, …) that make up the
Z-spectrum (normalised intensity ``I(\omega_\text{sat})/I_0``), plus the fitted
`R1`/`R2` rates and uncertainties. See
[Peak Lists and Output Files](peaklistformats.md) for the full format.

!!! note "Under development"
    `cest2d` is currently intended for data exploration and visualisation of Z-spectra.
    Peaks are fitted to a model of no exchange only (single-site relaxation); fitting to
    full exchange models (two-state Bloch-McConnell equations) is not yet implemented in
    the GUI and should be done in a separate analysis step.
