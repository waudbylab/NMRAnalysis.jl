# Heteronuclear NOE

The heteronuclear NOE reports on fast internal motion: rigid residues have an ¹H–¹⁵N NOE
near 0.8 at high field, and flexible ones lower or negative values. `hetnoe2d` measures it
as the ratio of peak amplitudes with and without ¹H saturation,

```math
\text{hetNOE} = \frac{I_\text{sat}}{I_\text{ref}}.
```

## Usage

Give the reference spectrum, then the saturated one:

```julia
using NMRAnalysis
hetnoe2d("11/pdata/231", "11/pdata/232")
```

Give lists of each to average over several pairs. Given all the planes in one list,
`hetnoe2d` reads which are saturated from the `hetnoe.saturated` annotation, or asks:

```julia
hetnoe2d(10:13)
hetnoe2d(10:13; saturation=[false, true, false, true])
```

### Example

The hetNOE of the calreticulin C-terminal domain:

```julia
reference, saturated = example("hetnoe2d_calr")
hetnoe2d(reference, saturated)
```

```@setup hetnoe2d
using NMRAnalysis
NMRAnalysis.screenshot("hetnoe2d.png") do
    hetnoe2d(example("hetnoe2d_calr")...)
end
```

![The hetnoe2d window, before any peaks are picked](hetnoe2d.png)

## Output

| Parameter | Meaning |
|---|---|
| `hetnoe` | I_sat / I_ref |
| `amp` | Amplitude in the reference spectrum |
