# Recovery Analysis

```@docs
recovery2d
```

The `recovery2d` function measures longitudinal relaxation from an inversion recovery or
saturation recovery experiment. Peak amplitudes are fitted to a magnetisation recovery
model:

```math
I(\tau) = A \left(1 - C \exp\!\left(-R\tau\right)\right)
```

where ``R`` is the recovery rate (s⁻¹), ``A`` is the equilibrium amplitude, and ``C``
is the recovery amplitude factor. For an ideal inversion recovery experiment ``C = 2``,
while for saturation recovery ``C = 1``.

![Screenshot from analysis with a magnetisation recovery model](../../assets/recovery2d.png)

## Usage

```julia
using NMRAnalysis

# Variable delay list
t = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
recovery2d("33/pdata/1", t)

# Read delays from a file
recovery2d("33/pdata/1", "vdlist.txt")
```

When the experiment is stored as a pseudo-3D dataset in a single processed directory, a
single path string is passed; when each delay is a separate experiment, pass a vector of
paths alongside the corresponding delay values.

## Output

Clicking **Save to folder** writes fitted parameters to `fit-results.txt`. Columns
include:

| Column | Description |
|--------|-------------|
| `label` | Peak label |
| `R_value` | Fitted recovery rate R (s⁻¹) |
| `R_uncertainty` | Uncertainty in R (s⁻¹) |
| `A_value` | Equilibrium amplitude A |
| `A_uncertainty` | Uncertainty in A |
| `C_value` | Recovery factor C |
| `C_uncertainty` | Uncertainty in C |
