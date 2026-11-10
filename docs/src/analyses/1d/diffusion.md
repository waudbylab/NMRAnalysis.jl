# Diffusion

`diffusion1d` fits the integral of each region against gradient strength with the
Stejskal-Tanner equation,

```math
I(g) = A \exp \left( -\left[ \gamma\delta\sigma g G_\mathrm{max} \right]^2
\left[\Delta - \delta/3 \right] D \right),
```

where ``g`` is the relative gradient strength, ``\delta`` the gradient pulse length,
``\Delta`` the diffusion delay, ``\sigma`` the gradient shape factor and ``\gamma`` the
gyromagnetic ratio. The hydrodynamic radius then follows from the Stokes-Einstein relation,

```math
D = \frac{k_B T}{6\pi \eta r_H},
```

with the solvent viscosity ``\eta`` at the sample temperature from Cho et al., *J. Phys.
Chem. B* **103**, 1991 (1999).

## Usage

Bruker doesn't store the gradient ramp, so give it as relative strengths from 0 to 1, one
per spectrum:

```julia
using NMRAnalysis
diffusion1d("106"; gradients=0.05:0.1:0.95)
```

Without `gradients`, you're asked for the first and last strengths, in percent, of a linear
ramp: the values you typed into TopSpin's `dosy` dialogue.

The other parameters are read from the acquisition parameters and offered for you to
confirm or change:

```
Gradient pulse length δ = 4000.0 µs (2*p30). Press enter to accept, or type a new value:
Diffusion delay Δ = 0.1 s (d20). Press enter to accept, or type a new value:
```

| Keyword | Meaning | Default |
|---|---|---|
| `gradients` | Relative gradient strengths, 0 to 1 | asked for |
| `δ` | Gradient pulse length, s | `2 × p30` |
| `Δ` | Diffusion delay, s | `d20` |
| `σ` | Gradient shape factor | 0.9 for `SMSQ`, 0.6366 for `SINE` |
| `Gmax` | Maximum gradient strength, T m⁻¹ | 0.55 |

The solvent and temperature come from the acquisition parameters. The hydrodynamic radius
is reported for `D2O` and `H2O+D2O` only.

### Example

An XSTE experiment on α-synuclein:

```julia
diffusion1d(example("diffusion1d_asyn"); gradients=0.05:0.1:0.95)
```

```@setup diffusion1d
using NMRAnalysis
res = NMRAnalysis.screenshot("diffusion1d.png") do
    reopen(exampleanalysis("diffusion1d_asyn"))
end
```

![The diffusion1d window, with the Stejskal-Tanner fit](diffusion1d.png)

## Output

| Parameter | Meaning |
|---|---|
| `A` | Amplitude |
| `D` | Diffusion coefficient, 10⁻¹⁰ m² s⁻¹ |
| `viscosity` | Solvent viscosity, mPa s |
| `rH` | Hydrodynamic radius, Å |

The example above returns:

```@example diffusion1d
res # hide
```
