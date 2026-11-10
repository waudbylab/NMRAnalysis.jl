# Titrations

`titration2d` fits a 1:1 binding isotherm to the chemical shift perturbations in a
titration series, with one dissociation constant shared by every residue, and the free and
bound shifts of each residue in each dimension.

For ``P + L \rightleftharpoons PL``, the observed shift moves from the free shift towards
the bound one in proportion to the fraction of protein bound, ``f``:

```math
\Delta\delta = f \,(\delta_\text{bound} - \delta_\text{free}).
```

Given the ligand concentrations ``L_\text{T}`` alone, the free ligand is taken to be the
total, and

```math
f = \frac{L_\text{T}}{K_\text{d} + L_\text{T}}.
```

Given the protein concentrations ``P_\text{T}`` as well, the exact solution is used, which
accounts for ligand depletion and dilution:

```math
f = \frac{1}{2P_\text{T}}\left(P_\text{T} + L_\text{T} + K_\text{d}
       - \sqrt{(P_\text{T} + L_\text{T} + K_\text{d})^2 - 4\,P_\text{T} L_\text{T}}\right).
```

## Usage

```julia
using NMRAnalysis
titration2d(1:11)
titration2d(1:11; L0=[0, 10, 20, 40, 80, 160], P0=[50, 49, 48, 46, 43, 38])
```

`L0` and `P0` are the total ligand and protein concentrations in each plane, in any unit, as
long as it's the same for both. ``K_\text{d}`` comes out in that unit. Without them, they're
read from the sample metadata of each experiment: the plane whose sample holds only one
molecule is the free protein, and the other molecule named is the ligand. The
concentrations found are printed as a table before the window opens. If the roles can't be
worked out, you're asked.

Peaks move from plane to plane, so they're added as in [Peak tracking](peaktracking.md):
`T` to add a peak and follow it through the planes, or `L` to place it along a line. Label
peaks with residue numbers, such as `L23`, for the summary plot.

### Example

The FIR RRM1–RRM2 domains titrated with the Nbox peptide. The concentrations are in the
sample metadata:

```julia
titration2d(example("titration2d_fir_nbox"))
```

```@setup titration2d
using NMRAnalysis
expt = NMRAnalysis.screenshot("titration2d.png"; peak="X14") do
    reopen(exampleanalysis("titration2d_fir_nbox"))
end
```

![The titration2d window, with the shift changes of the selected peak and the fitted isotherm](titration2d.png)

## Output

The global ``K_\text{d}`` is shown in the window and saved in `global.csv`. For each peak:

| Parameter | Meaning |
|---|---|
| `Xfree`, `Yfree` | Free shifts, ppm |
| `Xbound`, `Ybound` | Bound shifts, ppm |
| `dX`, `dY` | Shift changes on binding, ppm |
| `CSP` | Combined shift perturbation, ppm |

For the example above:

```@example titration2d
Main.resultstable(results(expt)) # hide
```

Label the peaks with residue numbers, and `summaryplot(expt)` plots these against residue
number.

The combined perturbation weights the two dimensions (Williamson, *Prog. Nucl. Magn. Reson.
Spectrosc.* **73**, 1 (2013)):

```math
|\Delta\delta| = \sqrt{(w_x\,\Delta\delta_x)^2 + (w_y\,\Delta\delta_y)^2}.
```

The `weights` keyword sets ``(w_x, w_y)``: `(1.0, 0.14)` by default, for ¹H and ¹⁵N, or
about `(1.0, 0.25)` for ¹H and ¹³C. They affect only this combined value, not the fit.
