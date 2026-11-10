# CPMG relaxation dispersion

CPMG (Carr-Purcell-Meiboom-Gill) relaxation dispersion experiments detect conformational
exchange processes on the microsecond-to-millisecond timescale. The effective transverse
relaxation rate R₂,eff depends on the CPMG pulse frequency νCPMG: residues undergoing
exchange show elevated R₂,eff at low νCPMG that decreases as νCPMG increases.

The effective relaxation rate is calculated from the intensity ratio:

```math
R_{2,\text{eff}}(\nu_\text{CPMG}) = -\frac{1}{T} \ln\!\frac{I(\nu_\text{CPMG})}{I_\text{ref}}
```

where ``T`` is the total CPMG relaxation period and ``I_\text{ref}`` is the
reference plane intensity (zero CPMG pulses).

## Usage

The input is a single pseudo-3D dataset where the first plane is the reference spectrum
and subsequent planes are recorded at increasing CPMG frequencies.

```julia
using NMRAnalysis

# Specify CPMG frequencies directly (Hz); use 0 for the reference plane
cpmg2d("11"; T=0.04, vCPMG=[0, 25, 50, 75, 100, 200, 500])

# Alternatively, specify the number of CPMG cycles per plane (vCPMG = ncyc / T)
ncyc = [0, 1, 2, 3, 4, 8, 20]
cpmg2d("11"; T=0.04, ncyc=ncyc)

# Or read the cycle numbers from the vclist
cpmg2d("11"; T=0.04)
```

If neither `vCPMG` nor `ncyc` is given, the cycle numbers are read from the `vclist`, and
if `T` or the cycle numbers can't be found you are asked for them. `skipplanes` leaves
planes out of the fit, as long as at least one reference plane remains.

!!! note
    When using `ncyc`, `vCPMG` is calculated automatically as `ncyc / T`. The
    reference plane must have `ncyc = 0` (or `vCPMG = 0`).

## Output

`results.csv` holds the reference rate `R20` (s⁻¹) of each peak, and `series.csv` the
amplitude in each plane, from which the R₂,eff profile follows.

!!! note "Under development"
    `cpmg2d` is currently intended for data exploration and extraction of R₂,eff
    profiles. Peaks are reported with a reference rate R₂⁰ only — no exchange model
    (fast-exchange, Bloch-McConnell, etc.) is fitted in the GUI. Exchange model fitting
    should be performed in a separate analysis step using the exported dispersion profiles.
