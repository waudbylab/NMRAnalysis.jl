# Peak Lists and Output File Formats

All 2D analyses read a peak list as input and write their results to a single
file. This page describes both, and the recommended (but not enforced) way to
label peaks.

## Recommended labelling

Any labelling convention is accepted — the parser never rejects a label. However, the
following conventions let the software extract residue numbers and atom types
for sorting and summary plots:

| Group | Pattern | Examples |
|-------|---------|----------|
| Backbone amides | one-letter code + residue number | `G10`, `S11`, `I13` |
| Methyls / sidechains | … + atom name | `I13CD1`, `L26CD1`, `L26CD2`, `V70CG1`, `M98CE` |
| Nucleic-acid atoms | … + atom name | `A12C8`, `G5C1'` |

Peaks whose label does not begin with a standard residue code — including the
default `X1`, `X2`, … given to newly picked, unassigned peaks — are assigned a
**negative** residue number and are omitted from residue-number summary plots by
default.

## Input: peak lists

A peak list tells the software where the peaks are. The only information the
reader needs is, for each peak, a **label** and its two **chemical shifts**
(`x` = direct/F1, `y` = indirect/F2, in ppm). Everything else is derived from
the label or measured by fitting.

The simplest input is a plain three-column text file, one peak per line:

```
G10   8.40   121.0
S11   9.02   115.7
I13   7.65   119.3
```

Fields may be separated by spaces, tabs, or commas. Lines beginning with `#`
are ignored. No header is required.

You can also load a previously saved `results.csv` (see below) to resume work or
to seed a new analysis from existing positions — the reader takes the `label`,
`x` and `y` columns and ignores the rest.

!!! note "Only label, x and y are read"
    When a file is loaded, **only the label and the two chemical shifts are
    used**. Residue number, residue type, atom name, linewidths, amplitudes and
    derived parameters are all ignored on input (the residue number and atom are
    re-derived from the label). You never have to reproduce those columns to
    re-use a file as input.

## Output: `results.csv`

Clicking **Save to folder** writes a single file, `results.csv`, with one row
per peak. Experiment metadata is written as `#`-comment lines, followed by an
ordinary header row and the data:

```
# Analysis type: Heteronuclear NOE
# Filename: /path/to/data
# Number of peaks: 3
label,resnum,resname,atom,x,x_err,y,y_err,R2x,R2x_err,R2y,R2y_err,amp[1],amp[1]_err,amp[2],amp[2]_err,hetnoe,hetnoe_err
G10,10,G,,8.40,0.01,121.0,0.05,30.1,1.2,15.2,0.8,4.5e5,2e3,3.6e5,2e3,0.78,0.04
```

Because the header is a real (uncommented) row, the file opens directly in a
spreadsheet and is read by, e.g., `pandas.read_csv("results.csv", comment="#")`.

The columns are:

| Column | Meaning |
|--------|---------|
| `label` | Peak label (editable in the GUI) |
| `resnum` | Residue number derived from the label (negative for unassigned peaks) |
| `resname` | One-letter residue code derived from the label |
| `atom` | Atom name derived from the label (blank for backbone amides) |
| `x`, `y` | Fitted chemical shifts (ppm), each with an `_err` uncertainty |
| `R2x`, `R2y` | Fitted linewidths in the direct/indirect dimensions (s⁻¹) |
| `amp[i]` | Fitted amplitude in spectrum/plane *i* |
| derived | Experiment-specific results (e.g. `hetnoe`, `R20`, `PRE`, `eta`, `R`) |

**Per-plane quantities are indexed in square brackets** (`amp[1]`, `amp[2]`, …).
A single-spectrum analysis has just `amp[1]`. The same `name[i]` convention would
extend to positions (`x[i]`) if a future analysis allowed peak positions to vary
between planes.

Each value column is immediately followed by its uncertainty (`value`,
`value_err`). The derived experiment parameters appear last, with the primary
result first. Existing files are backed up with an `.old` extension before being
overwritten.

`results.csv` is both the results table (for plotting and downstream analysis,
e.g. with [`summaryplot`](summary.md)) and a valid input file (for
reloading peak positions).

## Plotting summaries

[`summaryplot`](summary.md) plots a fitted parameter against residue number, from a
live experiment or one or more saved `results.csv` files. See the
[Summary Plots](summary.md) page for full details and examples.
