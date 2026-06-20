# Heteronuclear NOE

```@docs
hetnoe2d
```

The heteronuclear NOE (hetNOE) experiment measures cross-relaxation between ¹H and a
heteronucleus (typically ¹⁵N or ¹³C). It is widely used as a sensitive reporter of
fast internal dynamics: rigid residues typically show hetNOE values close to the maximum
(~0.8 for ¹H–¹⁵N at high field), while flexible regions show reduced or even negative
values.

The NOE is measured as the ratio of peak intensities in a saturated and a reference
spectrum:

```math
\text{hetNOE} = \frac{I_\text{sat}}{I_\text{ref}}
```

<!-- screenshot: assets/hetnoe2d.png -->

## Usage

Provide paired reference and saturated spectra. `reference` and `saturated` can
each be a single filename or a list of filenames; when lists are provided, results
are averaged across all pairs.

```julia
using NMRAnalysis

# Single reference / saturated pair
hetnoe2d("expno1/pdata/231", "expno1/pdata/232")

# Multiple pairs (results are averaged across pairs)
hetnoe2d(
    ["expno1/pdata/231", "expno2/pdata/231"],  # references
    ["expno1/pdata/232", "expno2/pdata/232"],  # saturated
)
```

## Output

Clicking **Save to folder** writes all results to `results.csv`. The derived
columns are:

| Column | Description |
|--------|-------------|
| `hetnoe`, `hetnoe_err` | Measured hetNOE = I_sat / I_ref and uncertainty |
| `amp`, `amp_err` | Reference amplitude and uncertainty |

See [Peak Lists and Output Files](fileformats.md) for the full format. Plot the
hetNOE against residue number with `summaryplot(expt)`.
