# Summary Plots

`summaryplot` plots a fitted parameter against residue number, from a live
experiment or one or more saved `results.csv` files. It is an ordinary Makie
figure using whichever backend is active, so it is interactive under GLMakie
and can be saved with `save("summary.pdf", fig)` under CairoMakie.

## Basic usage

```julia
fig = summaryplot(expt)                                   # default parameter, current experiment
fig = summaryplot("run1/results.csv"; param=:R20)         # specific parameter from a file
fig = summaryplot("run1/"; param=:R20)                    # folder containing results.csv
fig = summaryplot(["a/results.csv", "b/results.csv"]; param=:PRE)  # stacked panels
```

## Plot style

- **Backbone/amide peaks** (labels such as `A10N`, `G23HN`) → scatter of value
  vs residue number with error bars.
- **Atom-typed peaks** (e.g. methyls `I13CD1`, `L26CD2`) → bar chart ordered by
  `(residue, atom)` with peak-label ticks, so stereospecific pairs do not overlap.
  The style is chosen automatically per panel.
- Unassigned peaks (default `X#` names) are omitted unless every peak is
  unassigned, or `include_unassigned=true` is passed.

## Stacked panels

Passing a vector of sources produces vertically stacked panels, one per source.
By default each panel uses its own default parameter (so a mix of experiment
types — relaxation and hetNOE, for example — each show their own result).
Pass `param=:R2` to use the same parameter for every panel, or a vector such as
`param=[:R2, :hetnoe]` to set them individually.

```julia
# Stacked panels, each with its own default
fig = summaryplot(["relax/", "noe/"])

# Same parameter across all panels
fig = summaryplot(["wt/", "mutant/"]; param=:R20, ylabel="R₂⁰ / s⁻¹")

# Per-panel parameters
fig = summaryplot([expt_relax, expt_noe]; param=[:R2, :hetnoe])
```

## API reference

```@docs
summaryplot
```
