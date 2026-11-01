# Peak tracking

`peaktrack1d` follows the chemical shifts of peaks across a series of 1D spectra: a pH or
temperature series, a titration, or anything else where signals move from one spectrum to
the next. It draws the spectra as a stack, and you pick each peak by painting along its
trajectory with the mouse. A peak can be a singlet or a multiplet of several lines, and one
peak serves as the chemical shift reference.

Unlike the other 1D routines, `peaktrack1d` measures positions rather than integrals, and
fits nothing: the shifts are the result.

<!-- TODO screenshot: docs/src/assets/peaktrack1d-window.png
     The peaktrack1d window on a pH series: a stack of about twenty spectra, the reference
     picked in black, two singlets and one doublet picked in colour, the doublet's bars and
     dashed centre path visible, and the status line above the plot. -->

## Run the analysis

Call `peaktrack1d` with your experiments. Name the planes with `coordinates`, so that the
output says what distinguishes them, and name the reference:

```julia
using NMRAnalysis

pH = [5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
p = peaktrack1d(11:16; coordinates="pH" => pH, reference="DSS")
```

The window stays open until you close it, and `peaktrack1d` then returns the analysis.

### Choose how the spectra are processed

By default, `peaktrack1d` processes the raw data itself, so it doesn't depend on how the
spectra were processed in TopSpin. It phases every experiment automatically, and a slider
under the plot sets the exponential line broadening, which you can narrow to separate
overlapping lines while you pick.

To use the spectra processed in TopSpin instead, pass `process=false`:

```julia
p = peaktrack1d(11:16; process=false)
```

You can also pass spectra you have already loaded with NMRTools, either processed or as
FIDs from `loadfid`.

When processing, `phase` controls the phase correction:

| `phase` | Phases used |
|---|---|
| `:auto` (default) | found automatically for each experiment |
| `:pdata` | PHC0 and PHC1 from each experiment's processed data |
| a vector of `(p0, p1)` pairs | the ones you give, in degrees |
| a path to `processing.csv` | the ones saved by an earlier analysis |

In a series of three or more experiments, `peaktrack1d` checks the first-order phases
against each other. First-order phase comes from the acquisition rather than the sample, so
it should barely change along a series. A first-order phase far from the series median
usually means one line was phased correctly and another wrongly, so `peaktrack1d` sets it to
the median and fits the zero-order phase again. The REPL reports each experiment it
corrects, and `processing.csv` records it as `refit`.

## Pick a peak

The status line above the plot lists the peaks. The current peak is in brackets; the
reference comes first, and a `+` at the end is the vacant slot for a new peak.

1. Select a peak with `←` and `→`, or with the number keys, where `0` is the reference.
   To start a new peak, move to the `+` at the end.
2. Drag along the peak's trajectory through the stack. In each spectrum the brush crosses,
   the tallest peak inside the shaded area is picked.
3. Correct individual spectra: click to pick the spectrum under the cursor, drag a marker
   to move it, or right-click a marker to delete it. Shift-drag erases picks under the
   stroke.
4. Press `r`, type a name and press enter. A name can't contain an underscore or a comma: an
   underscore becomes a hyphen, and a comma is dropped.

If a strong neighbour is picked instead of your peak, narrow the brush with `[`. If two
lines overlap so that no maximum separates them, press `m` to turn snapping off and place
the picks by hand.

## Pick a multiplet

Pick each line of the multiplet as a separate peak, then group them:

1. Pick the first line, and name the peak.
2. Move to the `+` slot and pick the second line.
3. Press `g`, then click any marker of the first line. The second line joins the first
   line's peak, and keeps its name.

A grouped multiplet shows a bar across its lines in each spectrum, with a tick at its
centre and a dashed line joining the centres. Its chemical shift in each spectrum is the
mean of its lines, which is the centre of any symmetric first-order multiplet. A spectrum
in which one of the lines is missing has no centre, rather than a biased one.

Press `u` to move the current line back out into a peak of its own, or `d` to delete it.
The left and right keys step through every line of a multiplet in turn.

## Reference the shifts

Pick the reference like any other peak. Once it is picked in every spectrum, press `a` to
align the spectra so that the reference sits at `refshift` (0 ppm by default). Alignment
changes only the display. The shifts are always stored as measured, and the output gives
each centre referenced as well.

While the spectra are aligned, the reference can't be edited; press `a` again first.

## Save and reload

Press **Save** to write the output folder named in the box beside it, `out` by default.
Press **Load** to restore the picks saved there. To carry on from a saved analysis in a
new session, or on reprocessed data:

```julia
p = peaktrack1d(11:16; peaklist="out/peaklist.csv")
```

To reopen the window on an analysis you still have, pass it back:

```julia
p = peaktrack1d(p)
```

## Use the results

Show `p` in the REPL for a table of the centre of every peak in every spectrum. To paste
that table into a spreadsheet, press `c` in the window or run:

```julia
clipboard(p)
```

`planeresults(p)` gives every line in every spectrum as a vector of named tuples, the same
content as `series.csv`.

Press `s` in the window to save a publication figure of the region currently shown, at a
single journal column's width. For more control, draw it with `stackplot` and save it with
CairoMakie:

```julia
using CairoMakie
save("stack.pdf", stackplot(p; size=(500, 340)); backend=CairoMakie, pt_per_unit=1)
```

## Output files

| File | Contents |
|---|---|
| `summary.txt` | the processing, the reference, each peak's range of shifts, and the call that repeats the analysis |
| `peaklist.csv` | every pick as placed: what **Load** and `peaklist` read |
| `series.csv` | every line in every spectrum: its shift, height and signal-to-noise ratio, and its peak's centre and referenced shift |
| `shifts.csv` | the shift of every line as measured, one row per spectrum and one column per line |
| `referencedshifts.csv` | the same, referenced to the reference at `refshift`; written once the reference is picked in any spectrum |
| `heights.csv` | the height of every line, in the same layout |
| `processing.csv` | the phases, line broadening and zero-filling, when `peaktrack1d` processed the data |
| `stack.pdf` | the region shown in the window |
| `peaks/` | one CSV and one plot of shift against the first coordinate for each peak |

In `series.csv`, each line is identified by its peak's `label`, its `line` number
(counting from the highest shift) and the `plane`. The `centre (ppm)` and
`referenced (ppm)` columns describe the whole peak, so they repeat down its lines.

The three wide files are the ones to open in a spreadsheet or plot from. Their first columns
are `plane` and your coordinates. A singlet's column is named by its label, and the lines
of a multiplet by its label and line number, `MOPS[1]` and `MOPS[2]`, counting from the
highest shift. A blank pick is `NA`.

A height is the intensity at the line's maximum, divided by the number of scans and the
receiver gain, and by the sample concentration where one is recorded (NMRTools' `scale`),
so that spectra recorded with different settings compare. It is a height, not an integral or an amplitude: when a line broadens, as it often
does near a pKa or the midpoint of a titration, its height falls even though nothing has
been lost.

## Keys

| Key | Action |
|---|---|
| `←` `→` | previous or next line |
| `0` to `9` | select a peak; `0` is the reference |
| `r` | rename the current peak |
| `g`, then click | join the clicked line's peak to the current one |
| `u` | move the current line into a peak of its own |
| `d` | delete the current line |
| `a` | align the spectra to the reference |
| `m` | snapping on or off; shift inverts it while held |
| `↑` `↓` | vertical scale |
| `[` `]` | brush and search radius |
| `backspace` | undo |
| `p` | print the table in the REPL |
| `c` | copy the table |
| `s` | save a figure of the region shown |
| `q` | quit |

Scroll to zoom the chemical shift axis, right-drag to pan, and ctrl-click to reset the view.
