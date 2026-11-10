# Cross-correlated relaxation (CCR)

Cross-correlated relaxation between two mechanisms, such as the ¹H–¹⁵N dipolar coupling and
the ¹⁵N CSA, reports on bond orientation and dynamics. `ccr2d` measures the CCR rate η from
spectra in which a coherence decays or builds up over a relaxation period ``T``:

```math
\tanh(\eta T) = \frac{I_\text{buildup}}{I_\text{decay}}.
```

With symmetric reconversion, two buildup and two decay spectra cancel the effect of
auto-relaxation:

```math
\tanh(\eta T) = \sqrt{\frac{I_{\text{bu},1} \cdot I_{\text{bu},2}}{I_{\text{dec},1} \cdot I_{\text{dec},2}}}.
```

## Usage

Give the decay spectra, then the buildup, and the relaxation period `T` in seconds:

```julia
using NMRAnalysis
ccr2d(1, 2; T=0.08)                # one decay, one buildup
ccr2d([1, 4], [2, 3]; T=0.08)      # symmetric reconversion
```

Without `T`, you're asked for it.

### Example

Symmetric reconversion on the calreticulin C-terminal domain at 1.2 GHz, with
``T = 40`` ms:

```julia
decay, buildup = example("ccr2d_calr")
ccr2d(decay, buildup; T=0.04)
```

```@setup ccr2d
using NMRAnalysis
expt = NMRAnalysis.screenshot("ccr2d.png") do
    reopen(exampleanalysis("ccr2d_calr"))
end
Main.savepng("ccr2d-summary.png", summaryplot(expt))
```

![The ccr2d window](ccr2d.png)

## Output

| Parameter | Meaning |
|---|---|
| `eta` | CCR rate η, s⁻¹ |
| `amp` | Reference amplitude |

For the example above:

```@example ccr2d
Main.resultstable(results(expt)) # hide
```

![η against residue number](ccr2d-summary.png)
