# Relaxation (R₁, R₂)

`relaxation2d` fits the amplitude of each peak against the relaxation delay as an
exponential decay,

```math
I(T) = A \exp(-RT).
```

Whether ``R`` is R₁ or R₂ depends on the experiment. For an inversion or saturation
recovery, use [`recovery2d`](magnetisationrecovery.md).

## Usage

```julia
using NMRAnalysis
relaxation2d("20")                                   # one pseudo-3D experiment
relaxation2d(11:15; T=[0.01, 0.03, 0.06, 0.1, 0.2])    # one experiment per delay
```

The delays are found in this order:

1. `T`, in seconds, or a file holding one per line
2. `ncyc` × `cycletime`, for a delay counted in loops of `cycletime` seconds
3. the `relaxation.duration` annotation
4. the `vdlist`
5. the `vclist` × `cycletime`

Otherwise, you're asked for them.

### Example

¹⁵N R₂ of α-synuclein, measured with a CPMG train whose loop lasts 28.09 ms:

```julia
relaxation2d(example("relaxation2d_asyn"); cycletime=0.02809)
```

```@setup relaxation2d
using NMRAnalysis
expt = NMRAnalysis.screenshot("relaxation2d.png"; peak="V40") do
    reopen(exampleanalysis("relaxation2d_asyn"))
end
Main.savepng("relaxation2d-summary.png", summaryplot(expt))
```

![The relaxation2d window, with the decay of V40](relaxation2d.png)

![Picking and fitting peaks in a relaxation series](../../assets/relaxation-demo.mov)

## Output

**Save** writes these parameters for each peak to `results.csv`, the summary plot below to
`summary.pdf`, and a plot of each peak's fit to `peaks/`. See [Output](overview.md#Output)
for the other files.

| Parameter | Meaning |
|---|---|
| `R` | Relaxation rate, s⁻¹ |
| `A` | Amplitude |

For the example:

```@example relaxation2d
Main.resultstable(results(expt)) # hide
```

![R₂ against residue number](relaxation2d-summary.png)

Give the summary plot a label for your experiment with `ylabel`:

```julia
summaryplot("out"; ylabel="R₂ / s⁻¹")
```
