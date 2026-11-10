# Methyl CCR (S²τc)

`methylccr2d` measures the order parameter of the methyl symmetry axis, as ``S^2\tau_c``,
from ¹H–¹H cross-correlated relaxation in a methyl group (Sun, Kay and Tugarinov, *J. Phys.
Chem. B* **115**, 14878 (2011)). A buildup series ``I_a`` and a decay series ``I_b`` are
recorded at a series of delays ``T``, and their ratio is fitted to

```math
\left|\frac{I_a}{I_b}\right| =
\frac{C\,\eta\,\tanh\!\left(\sqrt{\eta^2+\delta^2}\,T\right)}
{\sqrt{\eta^2+\delta^2} - \delta\,\tanh\!\left(\sqrt{\eta^2+\delta^2}\,T\right)},
```

where ``\eta`` is the CCR rate, ``\delta < 0`` couples the fast- and slow-relaxing
coherences, and ``C`` is 3/4 for the triple-quantum experiment and 1/2 for the
double-quantum one. For an ideal methyl geometry with ``r_{HH} = 1.813`` Å,

```math
S^2\tau_c = \frac{\eta}{K}, \qquad
K = \frac{9}{40}\left(\frac{\mu_0}{4\pi}\right)^2 \frac{\gamma_H^4\,\hbar^2}{r_{HH}^6}
\approx 3.61\times10^{9}\ \text{s}^{-2},
```

so ``S^2\tau_c`` in ns is about ``0.277\,\eta``, with ``\eta`` in s⁻¹.

## Usage

Give the decay series, then the buildup:

```julia
using NMRAnalysis
methylccr2d(12, 11)
methylccr2d(12, 11; T=[0.002, 0.006, 0.012, 0.02, 0.03], C=1/2)   # double quantum
```

The delays are read from the buildup series' annotations or `vdlist`, or asked for. Give
them in seconds with `T`.

### Example

The triple-quantum experiment on ILV-labelled FLN5:

```julia
decay, buildup = example("methylccr2d_fln5")
methylccr2d(decay, buildup)
```

```@setup methylccr2d
using NMRAnalysis
expt = NMRAnalysis.screenshot("methylccr2d.png"; zoom=2, contour=-1) do
    reopen(exampleanalysis("methylccr2d_fln5"))
end
Main.savepng("methylccr2d-summary.png", summaryplot(expt))
```

![The methylccr2d window, with the buildup to decay ratio of the selected methyl](methylccr2d.png)

## Output

**Save** writes these parameters for each peak to `results.csv`, the summary plot below to
`summary.pdf`, and a plot of each peak's fit to `peaks/`. See [Output](overview.md#Output)
for the other files.

| Parameter | Meaning |
|---|---|
| `S2tc` | ``S^2\tau_c``, ns |
| `eta` | CCR rate η, s⁻¹ |
| `delta` | Coupling term δ, s⁻¹ |

For the example:

```@example methylccr2d
Main.resultstable(results(expt)) # hide
```

![S²τc for each methyl](methylccr2d-summary.png)
