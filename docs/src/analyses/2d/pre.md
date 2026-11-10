# Paramagnetic relaxation enhancement (PRE)

!!! warning "Under development"
    PRE fitting is under active development. Validate results carefully before you use
    them in a publication.

A paramagnetic centre, such as a spin label, a metal tag or a gadolinium cosolvent,
speeds up the relaxation of nearby nuclei. Use `pre2d` to measure the ¹H transverse
enhancement Γ₂ from spectra with and without the paramagnetic agent, and `lpre2d` to
measure the longitudinal enhancement Γ₁ as well, from spectra at a series of recycle
delays.

Both fit every peak across all of its spectra at once. Each peak has one diamagnetic
amplitude and diamagnetic linewidths, and the enhancement then sets its linewidths and
height in every spectrum:

- Γ₂ broadens the ¹H linewidth, and in an HMQC the multiple-quantum linewidth in the
  indirect dimension too.
- Γ₂ also lowers the peak height by ``\exp(-\Gamma_2 T)`` during the
  coherence transfer delays.
- For `lpre2d`, the height recovers between scans as
  ``1 - \exp(-(R_1 + \Gamma_1) \tau)``, where ``\tau`` is the recovery time.

For a single paramagnetic spectrum this gives the familiar intensity ratio
(Battiste & Wagner, *Biochemistry* **39**, 5355 (2000)):

```math
\frac{I_\text{para}}{I_\text{dia}} = \frac{R_2 \exp(-\Gamma_2 T)}{R_2 + \Gamma_2}
```

## Measure Γ₂ from a diamagnetic and a paramagnetic sample

Give the two spectra and the sequence type. Γ₂ is reported in s⁻¹.

```julia
using NMRAnalysis

pre2d("diamagnetic/pdata/1", "paramagnetic/pdata/1"; type=:hsqc)

# experiment numbers work too
pre2d(11, 12; type=:hmqc)
```

## Measure Γ₂ from a concentration series

For a paramagnetic cosolvent, give the spectra and the concentration of the agent in each,
with `0` for the diamagnetic reference. Γ₂ is then reported per unit of concentration, in
whatever unit you give, for example s⁻¹ mM⁻¹.

```julia
pre2d([10, 11]; concentrations=[0, 2], type=:hmqc)            # one Gd concentration
pre2d([10, 11, 12, 13]; concentrations=[0, 1, 2, 5], type=:hsqc)
```

## Measure Γ₁ and Γ₂ from a series of recycle delays

With spectra at several recycle delays, `lpre2d` fits the diamagnetic R₁ and both
enhancements together, so Γ₂ is free of the bias that incomplete recovery gives at short
delays. Give the diamagnetic and paramagnetic spectra as two lists, or one list with its
concentrations.

```julia
# spin-labelled protein: nine delays without and nine with the label
lpre2d(10:2:26, 11:2:27; type=:hmqc)

# 2 mM Gd: the nine diamagnetic experiments first, then the nine with Gd
lpre2d([10:18; 103:111]; concentrations=[fill(0, 9); fill(2, 9)], type=:hmqc)
```

The recovery time of each spectrum is its d1 plus its acquisition time, read from the
acquisition parameters. To set the recovery times yourself, pass `recycledelays`: the
values are used exactly as given, with nothing added. For two lists, you can give one list
of delays that applies to both.

R₁ and Γ₁ describe the recovery in the sequence you used. A selective or BEST-type
sequence recovers faster than a non-selective one, so compare them only between
experiments recorded the same way.

## Parameters

| Parameter | Description |
|-----------|-------------|
| `type` | `:hsqc` or `:hmqc`. Required. In an HMQC, ¹H is transverse during t₁, so Γ₂ also broadens the indirect dimension. |
| `concentrations` | Concentration of the paramagnetic agent in each spectrum, `0` for the reference. |
| `T` | Time (s) for which ¹H is transverse outside acquisition and t₁. The default `:auto` takes 1/¹J: 10.9 ms for ¹⁵N (92 Hz) and 8.0 ms for ¹³C (125 Hz, methyl). |
| `recycledelays` | `lpre2d` only. Recovery time (s) of each spectrum, used as given. |
| `skipplanes` | Spectra to leave out of the fit. At least two concentrations must remain. |
| `peaklist` | A saved peak list to load as the window opens. |

The `:auto` transfer time assumes a plain HSQC or HMQC, where ¹H is transverse for two
delays of 1/(2J). In a sensitivity-enhanced HSQC the time is longer and differs between
coherence pathways (Iwahara, Tang & Clore, *J. Magn. Reson.* **184**, 185 (2007)), so give
`T` explicitly, or use a sequence without sensitivity enhancement.

## Check the fits

Each peak's plot shows its fit:

- For a pair, cross-sections through both spectra, the diamagnetic and paramagnetic in
  different colours, so you can see both the broadening and the loss of height.
- For a concentration series, the peak height relative to its diamagnetic amplitude,
  against concentration, with the fitted curve.
- For `lpre2d`, the peak height against recovery time, a series and fitted curve for each
  concentration.

## Output

Clicking **Save to folder** writes the usual files (see
[Peak Lists and Output Files](peaklistformats.md)). `results.csv` has these columns:

| Column | Description |
|--------|-------------|
| `PRE` (`pre2d`), `Gamma2` (`lpre2d`) | Transverse enhancement Γ₂ |
| `Gamma1`, `R1` | `lpre2d` only: longitudinal enhancement Γ₁ and diamagnetic R₁ |
| `R2x0`, `R2y0` | Diamagnetic linewidths in the direct and indirect dimensions (s⁻¹) |
| `Irel` | For a single paramagnetic spectrum, its peak height relative to the diamagnetic one |

Rates carry the unit s⁻¹ for a diamagnetic and paramagnetic pair, and no unit for a
concentration series, since the concentration unit isn't known. In `series.csv`, `amp` is
the height measured in each spectrum at the linewidths the fit gives it, and `amp_fit` the
height the model predicts.

`summary.pdf` stacks panels against residue number: the fitted enhancements, each peak's
height in every paramagnetic spectrum relative to the diamagnetic reference (for `lpre2d`,
at the same recycle delay), and the diamagnetic R₂ and, for `lpre2d`, R₁. To plot a single
parameter instead, pass it to `summaryplot`:

```julia
expt = lpre2d(10:2:26, 11:2:27; type=:hmqc)
summaryplot(expt; param=:Gamma1)
```
