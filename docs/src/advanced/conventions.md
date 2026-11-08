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
  analysis.toml        the call that produced it, which `reopen` repeats
  peaklist.csv         what the user picked: 2D peaks and where they were placed, 1D
                       regions and the noise position
  results.csv          one row per entity: everything reported for it
  series.csv           the measurements, one row per entity per plane
  global.csv           parameters fitted once across every entity (when there are any)
  covariance.csv       covariance matrix of global.csv's fitted parameters (when jointly fitted)
  correlation.csv      the same fit's correlation matrix
  correlation.pdf      correlation.csv as a heatmap
  <overview>.pdf       fit.pdf in 1D, summary.pdf in 2D
  calibration.pdf      1D: the calibration curve of a multi-power nutation analysis
  calibration/         the analysis of a B₁ calibration fitted for this one
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
is absent where nothing is reported per entity, as in a kinetics run. `covariance.csv`,
`correlation.csv` and `correlation.pdf` appear only where that global fit jointly optimises
more than one parameter, since a single parameter has no covariance structure to report.
`calibration.pdf` appears only where a nutation analysis measured more than one power
level, one power level being a point rather than a curve; an analysis adds a file of its
own this way through `saveextras!`. A `calibration/` folder appears where an analysis
fitted a B₁ calibration on its own behalf: it is that calibration's own output folder,
written when the analysis that used it is saved and not before, so an abandoned fit leaves
nothing behind.

Saving starts from an empty folder: an existing one is moved aside to `<name>_previous`,
replacing any earlier backup. That is what keeps a peak or region deleted since the last
save from leaving its plot and its data behind, looking like part of the current result,
and it means a mistyped folder name never destroys what was there.

Cluster plots have no CSV of their own: a cluster is a group of overlapping peaks rather
than an entity with its own parameters, so it stays at the top level under its existing
`cluster_LABEL.pdf` name.

## Input and output are separate files

`peaklist.csv` records what the *user* specified: where each peak was placed, where each
region sits, the fitting radii, the noise position. It has one name in every analysis,
whatever the entities are, and every entry point takes it the same way, as `peaklist`: the
file, or a folder holding one, such as an earlier output folder. Everything else in
the folder records what the *fit* produced.

They are different things with different lives. A peak list is picked once and then reused:
on a second dataset, for a different kind of analysis, handed to a colleague, or imported
from elsewhere. A fitted position belongs to one particular fit. Reading fitted positions
back as the next run's starting point conflates the two, and for a moving-peak experiment
there is no single fitted position to read.

So the input file is what **Load** reads, and it is written on every save. It holds one row
per entity, or one row per entity per plane where positions were placed plane by plane, with
a blank `plane` key meaning the row applies to every plane.

In 1D, a region's row gives its bounds, `lo (ppm)` and `hi (ppm)`. The noise marker
travels as a region named `noise`, as wide as the
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

1D peak tracking (`peaktrack1d`) measures positions rather than integrals, and its entity is
a multiplet made of one or more lines. Its rows are keyed by `label`, `line` (numbered from
the highest shift) and `plane`, with the line's `shift (ppm)`, `height` and `snr`. The
multiplet's own `centre (ppm)` and `referenced (ppm)` describe the multiplet in that plane
rather than the line, so they are repeated down its lines, as `fitstatus` is in 2D. Nothing
is fitted, so there is no `results.csv`; the picks go in `peaklist.csv` (keyed by `label`,
`reference`, `line` and `plane`), and the processing of time-domain data, phases included,
in `processing.csv`, which can be passed back as `phase` to repeat it.

Peak tracking also writes three **wide** files, `shifts.csv`, `referencedshifts.csv` and
`heights.csv`, for opening in a spreadsheet: one row per plane, keyed by `plane` and the
coordinates, and one value column per line of every multiplet, named by the multiplet's
label (`MES (ppm)`) or by its label and line number (`MOPS[1] (ppm)`). Since their headers
are built from labels the user typed, labels never contain an underscore: one typed or read
from a peak list becomes a hyphen. The wide files repeat what is in `series.csv`, which
remains the complete record.

2D adds `fitstatus`, how the lineshape fit of that peak ended: `converged`, `maxiter`
(stopped at the iteration limit), `bound` (a position or linewidth ended at its limit),
`timeout` (stopped at the time limit, the previous values kept) or `unfitted`. It is a key,
repeated down a peak's rows, and so that a file is never silently unfinished it is written
in `results.csv` too.

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

2D adds `resnum`, `resname` and `atom`, derived from the label, and `fitstatus` (see
`series.csv` above). The region or peak
*positions* are not here: where an entity sits is something the user chose, and it lives in
`peaklist.csv`.

## `global.csv`

Long, because these parameter sets are small and structurally unlike each other.

```
parameter,value,error,unit
Kd,12.4,0.8,uM
```

An exchange fit adds `initial`, `fixed` and `atbound` columns, since every parameter of a
joint fit is global and what it started at, whether it moved, and whether it converged onto
a bound rather than an interior optimum are all part of the result.

## `summary.txt`

The human-readable record, and the one file where numbers are rounded. It carries the
package version and the date, the input filenames and titles, the sample information, the
region or peak definitions with the noise position and integration width, the acquisition
parameters actually used, and the key results formatted with units. Paths are written in
full here and in the CSV headers, since a saved file has to say where its data were; what
is shown on screen is shortened instead, having a panel's width to live within.

Where a region was measured under several conditions, its results are written as a block
per condition, under a heading naming it, followed by whatever was derived across them:

```
peak1
-----
11.11 dB
  Nutation frequency     2033.0 ± 7.4 Hz
  90° pulse              123.0 ± 0.45 µs
  B₁ inhomogeneity       7.583 ± 0.79 %

31.11 dB
  Nutation frequency     204.9 ± 0.66 Hz
  90° pulse              1220.0 ± 3.9 µs
  B₁ inhomogeneity       9.731 ± 0.53 %

B₁ inhom. (smallest)     7.583 ± 0.79 %
Linearity                0.9965 ± 0.0021
```

A region with one series has nothing to group, and is written as the flat block it always
was. The GUI's results panel is built from the same text, so the two cannot drift apart.

It ends by saying how to repeat the analysis:

```
Reproduce:
    reopen("out")
```

## `analysis.toml` and `reopen`

Every analysis saves the call that produced it as `analysis.toml`, beside `summary.txt`,
and [`reopen`](@ref) reads it back. The call is recorded with every parameter as it was
resolved, and with the `peaklist.csv` saved beside it, so reopening opens the window with the
same data, parameters and peaks, and the fit repeats:

```toml
# The analysis that wrote this folder, with paths relative to it.
# reopen("<this folder>") in NMRAnalysis.jl repeats it.
nmranalysis = "0.5.0"
call = '''
relaxation2d("../11/pdata/1";
             relaxationtimes=[0.01, 0.03, 0.06, 0.1, 0.2, 0.4],
             peaklist="./peaklist.csv")
'''
```

Paths are written relative to the folder, as `./name` inside it and `../name` outside it, so
the folder can be moved, copied or shipped as example data with its experiments and still be
reopened. Any string argument that named a file or folder when the analysis was run is
written this way, and so is an experiment number given as data; reading the record resolves
every string starting `./` or `../` against the folder.

The record is data, never run as code. `reopen` parses the call and accepts only literal
values (numbers, strings, symbols, `nothing`, and vectors, tuples, named tuples and pairs of
them) and only a function this package exports. Keywords passed to `reopen` replace those
recorded: `reopen("out"; peaklist=nothing)` repeats the analysis without its peaks.

A call has no record where its data were passed as `NMRData` rather than as paths, since
there is then nothing to point at.

`AnalysisCall`, `analysiscall`, `writereproduce` and `reopen` live in `src/calls.jl`, shared
by every module.

## Entry points

The data is positional; everything else is a keyword. Where an analysis takes its data in
roles, the roles are the positions: `hetnoe2d(reference, saturated)`,
`rdc2d(isotropic, aligned)`, `ccr2d(decay, buildup)`. A coordinate list is a keyword even
where it is required (a saturation list, a gradient ramp, the times of a kinetics run),
because it may instead be resolved from the experiment, and an optional positional argument
that is sometimes inferred reads badly. Data may be given as paths or as Bruker experiment
numbers, which `asexptpath` turns into paths.

What the user picked is passed the same way everywhere, as `peaklist`: a `peaklist.csv`, a
folder holding one, or the entities themselves (`Region`s in 1D).

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
what the experiment was, and each has a keyword that answers it, so that a saved fit can be
reopened without one.

## Implementation status

| Module | `summary.txt` | `analysis.toml` | `results.csv` | `series.csv` | `global.csv` | per-entity files |
|---|---|---|---|---|---|---|
| Analysis1D | yes | yes | yes | yes | nothing global yet | yes |
| GUI2D | yes | yes | yes | yes | yes | yes |
| Exchange1D | yes | yes | yes | yes | yes | yes |
| PeakTrack1D | yes | yes | nothing fitted | yes | nothing global | yes |
| R1rho | old format | no | no | no | no | no |

Every 2D routine now resolves its parameters through the same chain as the 1D routines,
from the helpers in `src/prompts.jl`, and returns the analysis when its window closes
(`results(expt)` and `planeresults(expt)` give the two tables as value ± error).

Still outstanding:

- Analysis1D has no `global.csv` because no 1D analysis currently fits anything across
  regions. `postfitglobal!` exists and will need somewhere to put its results when one does.
- `analysis.toml` records the resolved arguments but not *where* each came from (an
  annotation, the `vdlist`, a question), which needs the resolution chain to report its own
  winner rather than only its result.
