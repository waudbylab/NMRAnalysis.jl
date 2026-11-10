# TRACT

TRACT (Lee et al., *J. Magn. Reson.* **178**, 72 (2006)) estimates the rotational
correlation time τc of a protein from the relaxation rates of the TROSY and anti-TROSY
components of the amide ¹⁵N doublet. Each is fitted as an exponential decay, and the
difference gives the cross-correlated relaxation rate,

```math
\eta_{xy} = \frac{R_\mathrm{anti} - R_\mathrm{TROSY}}{2}
= f \left[ 4 J(0) + 3 J(\omega_N) \right],
\qquad J(\omega) = \frac{2}{5} \frac{\tau_c}{1 + \omega^2\tau_c^2},
```

which is solved for τc. ``f`` holds the dipolar and CSA constants: an N–H bond length of
1.02 Å, a ¹⁵N CSA of 160 ppm, and an angle of 17° between the two tensors.

TRACT assumes a rigid molecule that tumbles isotropically. Internal motion lowers η, so τc
is an effective value, and a lower bound for a flexible protein.

## Usage

Give the TROSY experiment, then the anti-TROSY:

```julia
using NMRAnalysis
tract("12", "13")
```

Both share one integration region, which starts on the amide envelope from 7.5 to
9.5 ppm. The delays are read from each experiment's `vdlist`, or give them with `T`:

```julia
tract("12", "13"; T=[0.0, 0.01, 0.02, 0.04, 0.08])
```

### Example

TRACT on calreticulin:

```julia
trosy, antitrosy = example("tract_calr")
tract(trosy, antitrosy)
```

```@setup tract
using NMRAnalysis
res = NMRAnalysis.screenshot("tract.png") do
    tract(example("tract_calr")...)
end
```

![The TRACT window, with the TROSY and anti-TROSY decays](tract.png)

## Output

| Parameter | Meaning |
|---|---|
| `A`, `R` | Amplitude and relaxation rate (s⁻¹) of each component |
| `etaxy` | Cross-correlated relaxation rate, s⁻¹ |
| `tauc` | Rotational correlation time, ns |

`etaxy` and `tauc` are reported with the TROSY result. The example above returns:

```@example tract
res # hide
```
