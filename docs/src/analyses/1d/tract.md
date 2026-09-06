# TRACT

TRACT (TROSY for Rotational Correlation Times) estimates the rotational correlation time
τc of a protein from the difference in relaxation rates between the TROSY and anti-TROSY
components of the ¹⁵N-¹H doublet. It needs two experiments, recorded as a pair, and gives
a global tumbling time in a few minutes of measurement rather than a full relaxation
analysis.

## Running the analysis

```julia
using NMRAnalysis

results = tract("12", "13")     # TROSY experiment, then anti-TROSY
```

Call `tract()` with no arguments to be asked for the two folders.

Both experiments are loaded into one dataset and share a single integration region, so
that the two rates are measured over exactly the same signal. The
[analysis window](overview.md#Using-the-analysis-window) opens with the region set to the
7.5 to 9.5 ppm amide envelope, which is what TRACT is normally measured over. Adjust it if
your spectrum calls for it, put the noise marker somewhere empty, and press **Save**.

The TROSY and anti-TROSY decays are plotted in separate colours, and the results panel
gives the two rates, the cross-correlated relaxation rate η, and τc.

<!-- TODO screenshot (4): docs/src/assets/tract1d-window.png
     tract on a TROSY/anti-TROSY pair: the amide envelope region, both decays fitted in
     their own colours below, and the panel showing the two rates with η and τc under a
     "TRACT results" heading. Replaces tract-regions.png and tract-fit.png. -->

## Relaxation delays

The delays are read from each experiment's own `vdlist`. Pass `tau` to give them
explicitly, for both experiments at once, and you are asked for them if neither is
available.

```julia
tract("12", "13"; tau=[0.0, 0.01, 0.02, 0.04, 0.08])
```

## What is fitted

An exponential decay to each of the two series,

```math
I(t) = A \exp(-Rt),
```

from which the cross-correlated cross-relaxation rate follows as

```math
\eta_{xy} = \frac{R_\mathrm{anti} - R_\mathrm{TROSY}}{2},
```

and τc by inverting

```math
\eta_{xy} = f \left[ \frac{4}{5}\tau_c + \frac{3}{5}\frac{\tau_c}{1 + \omega_N^2\tau_c^2} \right],
```

where ``f`` collects the dipolar and CSA constants for the amide group: the N-H bond
length (1.02 Å), the ¹⁵N chemical shift anisotropy (160 ppm), the angle between the two
tensors (17°), and the static field.

| Parameter | Meaning |
|---|---|
| `A` | Amplitude of each decay |
| `R` | Relaxation rate of each component, s⁻¹ |
| `ηxy` | Cross-correlated relaxation rate, s⁻¹ |
| `τc` | Rotational correlation time, ns |

η and τc describe the pair rather than either component, so they are recorded once, on the
TROSY series:

```julia
results = tract("12", "13")
trosy = first(r for r in results if r.group.which == :trosy)
param(trosy, :τc)     # ns, with uncertainty
```

## Assumptions

!!! warning "τc is an effective value"
    TRACT assumes isotropic tumbling of a perfectly rigid molecule. Internal motion,
    flexible termini and disorder all reduce the measured η, so a reported τc is a lower
    bound for a molecule that is not rigid. Read it as an effective correlation time, not
    a hydrodynamic one.

TRACT experiments are not currently recognised by [`analyse`](../analyse.md), because they
need a pair of files rather than one. Call `tract` directly.

## Citation

Lee, D., Hilty, C., Wider, G., & Wüthrich, K. (2006). Effective rotational correlation
times of proteins from NMR relaxation interference. *Journal of Magnetic Resonance*,
**178**, 72-76. doi: 10.1016/j.jmr.2005.08.014
