# Cross-Correlated Relaxation (CCR)

Cross-correlated relaxation (CCR) experiments measure the interference between two
relaxation mechanisms — for example, between dipole-dipole coupling and chemical shift
anisotropy — and are used to determine bond vector orientations, order parameters, and
rotational correlation times.

The CCR rate η is extracted by comparing the buildup and decay of spin-state-selective
coherences:

```math
\tanh(\eta T) = \frac{I_\text{buildup}}{I_\text{decay}}
```

For symmetric reconversion experiments, where two buildup and two decay spectra are
recorded to suppress contributions from auto-relaxation:

```math
\tanh(\eta T) = \sqrt{\frac{I_{\text{bu},1} \cdot I_{\text{bu},2}}{I_{\text{dec},1} \cdot I_{\text{dec},2}}}
```

<!-- screenshot: assets/ccr2d.png -->

## Usage

### Single decay / buildup pair

```julia
using NMRAnalysis

ccr2d("decay/pdata/1", "buildup/pdata/1", 0.08)
```

### Symmetric reconversion (two pairs)

```julia
ccr2d(
    ["decay1/pdata/1", "decay2/pdata/1"],
    ["buildup1/pdata/1", "buildup2/pdata/1"],
    0.08
)
```

The third argument `T` is the relaxation time constant in seconds during which the CCR
rate acts.

## Output

Clicking **Save to folder** writes all results to `results.csv`. The derived
columns are:

| Column | Description |
|--------|-------------|
| `eta`, `eta_err` | Fitted CCR rate η (s⁻¹) and uncertainty |
| `amp`, `amp_err` | Reference amplitude and uncertainty |

See [Peak Lists and Output Files](peaklistformats.md) for the full format. Plot η
against residue number with `summaryplot(expt)`.

## References

- Reif, B., Hennig, M. & Griesinger, C. (1997) Direct measurement of angles between
  bond vectors in high-resolution NMR. *Science* **276**, 1230–1233.
  [doi:10.1126/science.276.5316.1230](https://doi.org/10.1126/science.276.5316.1230)

- Pelupessy, P., Chiarparin, E., Ghose, R. & Bodenhausen, G. (1999) Simultaneous
  determination of φ and ψ angles in proteins from measurements of cross-correlated
  relaxation effects. *J. Biomol. NMR* **13**, 375–380.
  [doi:10.1023/A:1008365915981](https://doi.org/10.1023/A:1008365915981)

- Richter, C., Griesinger, C., Felli, I., Cole, P. T., Varani, G. & Schwalbe, H.
  (1999) Determination of sugar conformation in large RNA oligonucleotides from
  analysis of dipole–dipole cross correlated relaxation by solution NMR spectroscopy.
  *J. Biomol. NMR* **15**, 241–250.
  [doi:10.1023/A:1008319916075](https://doi.org/10.1023/A:1008319916075)
