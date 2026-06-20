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

Provide paired reference and saturated spectra. Each reference spectrum must be
immediately followed by its saturated counterpart in the file list, and the corresponding
`saturationlist` entry must be `false` (reference) or `true` (saturated).

```julia
using NMRAnalysis

# Single reference / saturated pair
hetnoe2d(
    ["ref/pdata/1", "sat/pdata/1"],
    [false, true]
)

# Multiple pairs (results are averaged across pairs)
hetnoe2d([
    "expno1/pdata/231",   # reference
    "expno1/pdata/232",   # saturated
    "expno2/pdata/231",   # reference
    "expno2/pdata/232",   # saturated
], [false, true, false, true])
```

## Output

Clicking **Save to folder** writes hetNOE values to `fit-results.txt`:

| Column | Description |
|--------|-------------|
| `label` | Peak label |
| `hetNOE_value` | Measured hetNOE = I_sat / I_ref |
| `hetNOE_uncertainty` | Propagated uncertainty |
