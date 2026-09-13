# Output and Interface Conventions

Every analysis in NMRAnalysis.jl writes the same set of files, with the same column rules,
whatever kind of experiment it was. This page is the specification those files are written
against, and is aimed at anyone adding or changing an analysis; the per-experiment pages
describe what each one means. [How an Analysis Works](pipeline.md) describes the three
stages that produce them.

The point is that a results folder should be readable without knowing which routine
produced it, and that a script written against one analysis should work against another.

## The output folder

```
out/
  summary.txt          human-readable record of the whole analysis
  peaklist.csv         2D: the peaks the user picked, and where they were placed
  regionlist.csv       1D: the regions the user picked, and the noise position
  results.csv          one row per entity: everything reported for it
  series.csv           the measurements, one row per entity per plane
  global.csv           parameters fitted once across every entity (when there are any)
  <overview>.pdf       fit.pdf in 1D, summary.pdf in 2D
  regions/             1D: one file pair per region
    signal.csv
    signal.pdf
  peaks/               2D: one file pair per peak
    L23N.csv
    L23N.pdf
  cluster_*.pdf        2D only: one plot per group of overlapping peaks
```

A per-entity CSV and its plot share a basename and differ only in extension, so the data
behind any plot is beside it. The per-entity CSVs are the rows of `series.csv` filtered to
that entity; the duplication is deliberate, since opening one peak's data should not
require filtering a file of several thousand rows.

A file is written only when it has something to say. `global.csv` appears only where the
analysis fits something across every entity, which today means a titration `Kd` and an
exchange `kex`; a relaxation fit has nothing global and the file is absent. `results.csv`
is absent where nothing is reported per entity, as in a kinetics run.

Saving starts from an empty folder: an existing one is moved aside to `<name>_previous`,
replacing any earlier backup. That is what keeps a peak or region deleted since the last
save from leaving its plot and its data behind, looking like part of the current result,
and it means a mistyped folder name never destroys what was there.

Cluster plots have no CSV of their own: a cluster is a group of overlapping peaks rather
than an entity with its own parameters, so it stays at the top level under its existing
`cluster_LABEL.pdf` name.

## Input and output are separate files

`peaklist.csv` and `regionlist.csv` record what the *user* specified: where each peak was
placed, where each region sits, the fitting radii, the noise position. Everything else in
the folder records what the *fit* produced.

They are different things with different lives. A peak list is picked once and then reused:
on a second dataset, for a different kind of analysis, handed to a colleague, or imported
from elsewhere. A fitted position belongs to one particular fit. Reading fitted positions
back as the next run's starting point conflates the two, and for a moving-peak experiment
there is no single fitted position to read.

So the input file is what **Load** reads, and it is written on every save. It holds one row
per entity, or one row per entity per plane where positions were placed plane by plane, with
a blank `plane` key meaning the row applies to every plane.

The 1D noise marker travels in `regionlist.csv` as a region named `noise`, as wide as the
widest signal region. Only its centre is read back; the width used to estimate a region's
uncertainty always matches that region's own.

2D peak lists in **Sparky** format can be imported. Sparky carries one position per peak
and nothing else, so it is an import format rather than the one written here: a trajectory,
the fitting radii and the uncertainties have nowhere to go in it. Its `w1`/`w2` columns are
matched to the direct and indirect axes by where the shifts actually fall, since reading
them in order would transpose every peak in a list written the other way round.

## Entities and coordinates

An **entity** is what an analysis reports on: a peak in 2D, a region in 1D, the whole
problem in a joint exchange fit. Every entity has a `label`.

A **coordinate** is anything that distinguishes one plane from another: a relaxation delay,
a gradient strength, a spin-lock field, a concentration, or a categorical tag such as TROSY
versus anti-TROSY. Coordinates may be numeric or categorical, and there may be several.

Which file a number goes in follows from what it describes, not from the stage that
produced it. Anything varying plane by plane is in `series.csv`; anything describing one
entity is on that entity's single row of `results.csv`, whether the fit produced it or a
later step derived it; anything shared by every entity is in `global.csv`. TRACT's τc is
computed after the fit but describes a region, so it is on the region's row; a titration
`Kd` is one number for every peak, so it is global.

## Column rules

These hold in every CSV.

Column names are **ASCII**, so that `df.tauc` works: `etaxy`, `tauc`, `pulse90`, `R1rho`.
The typeset names belong in the GUI and in `summary.txt`, not in a file header.

A column carrying a physical quantity names its **unit in parentheses** after a space, in
ASCII: `R (s-1)`, `tauc (ns)`, `D (1e-10 m2/s)`, `pulse90 (us)`. A dimensionless quantity
has no parentheses.

A value column `X` may be accompanied by an **uncertainty** column `X_err` and a **fitted
value** column `X_fit`, each repeating the unit: `R (s-1)`, `R_err (s-1)`. Repeating it
keeps the two symmetrical for anything reading them.

Every other column is a **key**. So the rule for reading one of these files generically is:
strip the parenthesised unit from each header, then any column whose name is `X`, `X_err`
or `X_fit` for some `X` is a value, and everything else identifies the row.

A **blank key** means the row applies to every value of that key. A `peaklist.csv` row with
`plane` left blank places that peak at the same position in every plane.

A parameter name **carries no underscore of its own**, so an underscore in a header always
separates the quantity from something added to it: the `_err` and `_fit` markers, or the
series a fitted parameter came from. `R_trosy` is the quantity `R`, measured in the series
tagged `trosy`, and that is how its label and unit are found.

Numbers are written at **full precision**. These are machine files and rounding is
irreversible; `summary.txt` is where numbers are rounded for reading. A value that does not
exist for a row is written `NA`, which is distinct from a blank key.

## `series.csv`

One row per entity per plane, in long form, because an analysis may have several
coordinates and several datasets.

```
source,label,plane,which,time (s),I,I_err,I_fit
12/pdata/1,amide,1,trosy,0.000,9421.3,12.4,9430.1
12/pdata/1,amide,2,trosy,0.005,8684.0,12.4,8687.5
13/pdata/1,amide,5,anti,0.000,9388.2,12.4,9401.7
```

`source` names the dataset each row came from and is always present, even when it is
constant. It is provenance, not meaning: two files may be replicates with identical
coordinates, so whatever physical variable distinguishes datasets (a concentration, a
spin-lock field, a TROSY tag) gets its own coordinate column as well.

`plane` is always present for the same reason and is load-bearing where the planes share
one file: every row of a pseudo-2D experiment has the same `source`, and the plane index is
then the only thing telling two rows apart.

`I_fit` is the model evaluated at the measured coordinates, not on a fine grid, so that
residuals are a subtraction. It is `NA` where nothing was fitted. A smooth curve for a
figure is recoverable from the parameters.

## `results.csv`

One row per entity, wide, because this is the table you sort by residue number and plot
against it.

```
label,tauc (ns),tauc_err (ns),etaxy (s-1),etaxy_err (s-1),A_trosy,A_trosy_err,R_trosy (s-1),R_trosy_err (s-1),A_anti,A_anti_err,R_anti (s-1),R_anti_err (s-1)
amide,14.86,0.42,19.04,0.51,9421.3,12.4,16.02,0.41,9388.2,12.4,54.10,0.93
```

Where an entity was measured under more than one condition, each condition's fitted
parameters are named apart on that same row rather than splitting the entity across rows.
That is what lets a quantity combining conditions, like TRACT's τc, sit beside the two rates
it came from with nothing left blank. The experiment's headline parameter
(`primaryparam`) comes first.

2D adds `resnum`, `resname` and `atom`, derived from the label. The region or peak
*positions* are not here: where an entity sits is something the user chose, and it lives in
`regionlist.csv` or `peaklist.csv`.

## `global.csv`

Long, because these parameter sets are small and structurally unlike each other.

```
parameter,value,error,unit
Kd,12.4,0.8,uM
```

An exchange fit adds `initial` and `fixed` columns, since every parameter of a joint fit is
global and what it started at and whether it moved are part of the result.

## `summary.txt`

The human-readable record, and the one file where numbers are rounded. It carries the
package version and the date, the input filenames and titles, the sample information, the
region or peak definitions with the noise position and integration width, the acquisition
parameters actually used, and the key results formatted with units.

It should also carry the Julia call that would repeat the analysis, together with where
each resolved parameter came from, which makes the annotation lookup auditable:

```
Reproduce:
    relaxation1d("11";
                 tau=[0.01, 0.03, 0.06, 0.10, 0.20, 0.40],
                 model=:exponential,
                 integration=(peakppm=8.21, noiseppm=-1.00, ppmwidth=0.60))

tau from vdlist; model from annotation relaxation.model; region selected interactively.
```

## Entry points

The data is positional; everything else is a keyword. A coordinate list is a keyword even
where it is required, because it may instead be resolved from the experiment, and an
optional positional argument that is sometimes inferred reads badly.

Keyword names are chosen to be informative rather than to match annotation keys, so
`relaxationtimes`, `Trelax` and `Tsat` keep the names a spectroscopist would use.

Where a parameter can be resolved rather than typed, it is looked for in a fixed order:
an explicit argument, then a pulse-sequence annotation, then sample metadata, then a Bruker
acquisition parameter, and finally a question asked before any window opens. `prompt=false`
(the default outside an interactive session) turns that last step into either a documented
default or an error naming the argument to pass.

Not every analysis resolves everything. An experiment combining many datasets, each with
its own lists of offsets, powers and delays, is not something a user can reasonably type,
so `exchange1d` reads those from annotations and does not offer a fallback. Its questions
are about the fit (which model, which molecules, which parameters to fix) rather than about
what the experiment was.

## Implementation status

| Module | `summary.txt` | `results.csv` | `series.csv` | `global.csv` | per-entity files |
|---|---|---|---|---|---|
| Analysis1D | yes | yes | yes | nothing global yet | yes |
| GUI2D | yes | yes | yes | yes | yes |
| Exchange1D | yes | yes | yes | yes | yes |
| R1rho | old format | no | no | no | no |

Still outstanding:

- GUI2D has no `Reproduce:` block, and no equivalent of `AnalysisCall` to record what each
  routine was given.
- Analysis1D has no `global.csv` because no 1D analysis currently fits anything across
  regions. `postfitglobal!` exists and will need somewhere to put its results when one does.
- The `Reproduce:` line records the resolved arguments but not *where* each came from (an
  annotation, the `vdlist`, a question), which needs the resolution chain to report its own
  winner rather than only its result.
