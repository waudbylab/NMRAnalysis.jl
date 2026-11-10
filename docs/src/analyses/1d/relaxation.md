# Relaxation (R₁, R₂)

`relaxation1d` fits the integral of each region against the relaxation delay, as an
exponential decay,

```math
I(T) = A \exp(-RT),
```

or as an inversion or saturation recovery,

```math
I(T) = A \left[ 1 - C \exp(-RT) \right],
```

where ``C`` is 2 for a complete inversion and 1 for saturation.

## Usage

```julia
using NMRAnalysis
relaxation1d("11")
```

The delays are read from the `relaxation.duration` annotation or the `vdlist`, and the model
from the `relaxation.model` annotation. Anything not found is asked for. To give them
yourself:

```julia
relaxation1d("11"; T=[0.01, 0.03, 0.06, 0.1], model=:recovery)
```

| Keyword | Meaning |
|---|---|
| `T` | Relaxation delays in seconds, one per spectrum |
| `model` | `:exponential` or `:recovery` |
| `peaklist` | Regions to start from: a saved `peaklist.csv`, or the folder holding one |
| `integration` | A single region, `(peakppm, noiseppm, ppmwidth)`, to fit without the window |

For an inversion recovery, step to the last spectrum before you place the regions: the
signal is largest there.

### Example

A ¹H T₂ series on a mixture of amino acids:

```julia
reopen(exampleanalysis("relaxation1d_aas"))
```

```@setup relaxation1d
using NMRAnalysis
res = NMRAnalysis.screenshot("relaxation1d.png") do
    reopen(exampleanalysis("relaxation1d_aas"))
end
```

![The relaxation1d window, with a T₂ decay fitted](relaxation1d.png)

## Output

| Parameter | Meaning |
|---|---|
| `A` | Amplitude |
| `R` | Relaxation rate, s⁻¹ |
| `C` | Recovery factor (recovery model only) |

The example above returns:

```@example relaxation1d
res # hide
```
