# 2D Experiment Analysis

NMRAnalysis.jl provides interactive GUIs for analysing 2D NMR experiments, including
relaxation, exchange, and NOE measurements. All functions follow the same pattern: they
load one or more processed 2D spectra, open an interactive window for peak picking and
fitting, and export results to a folder of your choice. When you close the window, the
function returns the analysis, so you can carry on with it in Julia:

```julia
expt = relaxation2d("11/pdata/1")
r = results(expt)           # one row per peak, each parameter as value ± error
[row.R for row in r]
planeresults(expt)          # one row per peak per plane
```

Any parameter you don't give is read from the pulse-sequence annotations or the acquisition
parameters where possible, and otherwise you are asked for it before the window opens. Pass
`prompt=false` to raise an error instead, for example in a script.

Every routine takes `skipplanes`, a list of planes to leave out of the fitting. Skipped
planes are still loaded and displayed, and their amplitudes are still measured, at the
lineshapes the other planes determine.

The [experiment-specific pages](fit.md) describe the available functions
and the theory behind each analysis. This page covers the shared GUI features common to
all of them.

![Screenshot of peak tracking](../../assets/peaktrack-demo.mov)

![Screenshot of relaxation fitting](../../assets/relaxation-demo.mov)

## Adding and Managing Peaks

Peaks are picked interactively using the mouse and keyboard. Move the cursor over a peak
in the contour plot to work with it.

| Action | Key / Button |
|--------|--------------|
| Add peak at cursor position | `A` |
| Track peak from cursor position (moving-peak experiments) | `T` |
| Place a peak at the maximum along a line in every plane (moving-peak experiments) | Hold `L` at one end, drag to the other, and release; or press `L` at each end |
| Widen or narrow the selected peak's fitting radii | `Shift` + `←` / `→` (x), `Shift` + `↑` / `↓` (y) |
| Return the selected peak to the default radii | `=` |
| Continue fits that stopped at the time or iteration limit | `C` |
| Continue them, without a time limit, until they converge | `Shift` + `C` |
| Cancel the fit in progress | `Esc` |
| Delete the selected peak | `D` or **Delete peak** button |
| Rename the selected peak | `R` or **Rename peak** button |
| Navigate to previous spectrum slice | `←` or **←** button |
| Navigate to next spectrum slice | `→` or **→** button |
| Raise contour base level | `↑` or **contour ↑** button |
| Lower contour base level | `↓` or **contour ↓** button |
| Reset axis zoom | **reset zoom** button |
| Show or hide the fitted lineshape overlay | **Fitting** toggle |
| Show or hide other spectra (moving-peak experiments, e.g. titrations, RDCs) | `S` or **(S)how all** toggle |
| Open a summary plot of the current results | **Summary plot** button (enabled once peaks are present) |
| Load a previously saved peak list | **Load peak list** button, or `peaklist=` when you start (see below) |
| Save all results to a folder | **Save to folder** button |
| Close the GUI window | `Q` or **(Q)uit** button |

Peak lineshapes are fitted in real time as you add or move peaks. The right panel shows
cross-sections (or a model fit plot, for relaxation-type experiments) for the currently
selected peak.

## Visual Feedback

The window background changes colour to indicate the current interaction mode:

| Background | Mode |
|------------|------|
| White | Normal |
| Salmon / orange | Fitting in progress |
| Light yellow | Adding a peak plane by plane, or marking a line |
| Light blue | Renaming a peak |
| Pale green | Moving a peak |
| Grey | Saving results to folder |

Peak markers are colour-coded:

| Colour | Meaning |
|--------|---------|
| Blue | Fitted peak |
| Red | Peak waiting to be fitted |
| Orange | Peak whose last fit didn't converge (see below) |
| Green | Currently selected peak |

## How peaks are fitted

Overlapping peaks are fitted together, as a cluster. Each peak is a 2D lineshape with a
position and linewidth in each dimension, truncated at its fitting radius. Its amplitude
in each plane is solved for exactly at every step of the fit, so for a series of fixed
peaks only the positions and linewidths are optimised, however many planes there are.
With more than one Julia thread (start Julia with `julia -t auto`), clusters are fitted in
parallel. The label beside the **Fitting** toggle shows progress through the clusters.

In a series of fixed peaks, each amplitude's uncertainty is that of the linear
least-squares fit of the amplitudes, with the positions and linewidths held at their fitted
values. It depends on the peak's shape, on its overlap with neighbouring peaks and on the
residual of the fit, not on the spectrum's noise level alone. An error in those shared
shapes would scale every amplitude by the same factor, which cancels in a relaxation rate,
a heteronuclear NOE, a CCR ratio or a CEST profile, so it isn't added to each plane. It
does affect an absolute amplitude, or a fitted prefactor such as `A`, which can therefore
be slightly underestimated where a linewidth is poorly determined; a `bound` status on
the peak shows when that is likely.

A fit that doesn't reach a converged optimum is flagged, and its peaks turn orange:

- a fit that runs past 30 seconds is stopped, leaving the previous values in place;
- a fit that reaches the iteration limit stops where it got to;
- a fit that ends with a position or linewidth at its limit is reported as such.

The info panel says which applies to the selected peak, and the label beside the
**Fitting** toggle counts the unfinished peaks. Press `C` to continue the stopped fits from
where they got to, with five minutes each, or `Shift` + `C` to keep continuing them, with no
time limit, round after round, until they converge or a round no longer moves them. The
REPL reports each round. `Esc` cancels either.

Positions are bounded to within the peak's radius of where you placed it. Linewidths are
bounded between 1 s⁻¹ and the broadest line the fitting window can determine, one whose full
width at half height spans twice the window (four radii), or 100 s⁻¹ if that is larger. For
a 0.04 ppm ¹H radius at 600 MHz that is about 300 s⁻¹. To fit broader peaks, widen their
radius. A position at its limit usually means the peak
needs moving or its radius widening. The status of each peak is saved in the `fitstatus`
column of `results.csv` and `series.csv`, and `summary.txt` lists any unfinished peaks.

## Recommended Workflow

1. Launch the appropriate analysis function with your input files.
2. Navigate to a representative spectrum using `←` / `→` or the slice slider.
3. Adjust contour levels with `↑` / `↓` until peaks are clearly visible.
4. Add peaks with `A` at each resonance you want to track.
5. Optionally rename peaks with `R` to match residue assignments.
6. Navigate through all slices to verify fit quality across the series.
7. Click **Save to folder** to write all output files to a chosen directory.

!!! tip
    For large datasets, it is efficient to pick peaks on one representative slice first,
    then step through remaining slices to check that the fits are good.

## Output Files

Clicking **Save to folder** writes the following files:

| File | Contents |
|------|---------|
| `summary.txt` | The record to read: where the data came from, the fitting radii, and the headline parameter for every peak |
| `peaklist.csv` | The peaks you picked and where you placed them, with the fitting radii. This is what **Load** reads |
| `results.csv` | One row per peak: its identity and the derived parameters (relaxation rates, NOE values, …), each with an uncertainty. This is the table to plot against residue number |
| `series.csv` | The measurements, one row per peak per plane: the plane index and its own coordinate, then the amplitude, the fitted amplitude, the position and the linewidths |
| `global.csv` | Anything fitted once across every peak, such as a titration `Kd`. Absent when there is nothing global |
| `summary.pdf` | Summary plot of the primary fitted parameter against residue number (or atom for methyl/non-backbone experiments) |
| `peaks/LABEL.pdf` and `.csv` | Each peak's fit plot, and the data behind it under the same name |
| `cluster_LABEL.pdf` | Zoomed 2D contour plot (first plane) with fitted lineshapes for each group of overlapping peaks |

Anything that varies plane by plane is in `series.csv`: the amplitude, the position and
the linewidths. A peak that does not move simply repeats its position down the rows, so the
layout is the same whether or not the peaks track. Each row names the plane's own
coordinate, so a relaxation series records the delay and a titration the concentration,
rather than an `amp[7]` column whose meaning has to be remembered.

`amp_fit` is the model evaluated at that plane's coordinate, so a residual is a
subtraction. It is `NA` where nothing is fitted through the amplitudes themselves: a
heteronuclear NOE, a CCR rate and a CEST profile are fitted from ratios or from
transformed intensities.

Every CSV has experiment metadata in `#`-comment lines above an ordinary header
row, so it opens directly in spreadsheets and `pandas`. Column headers carry
units in ASCII (`R (s-1)`, `x (ppm)`), and an uncertainty column repeats the unit
of the value it belongs to. An existing output folder is moved aside to
`<name>_previous` before saving, so each save starts clean and a peak you deleted
since the last one does not leave its plot behind. See
[Peak Lists and Output Files](peaklistformats.md) for the full column description.

## Loading and Resuming Analysis

The **Load peak list** button restores peak positions and labels from a saved
`peaklist.csv`, a Sparky peak list, or a simple `label x y` text file, so you can resume
work later or seed a new analysis from existing positions. Where peaks were tracked plane
by plane, the whole trajectory is restored, as are any radii you set for individual peaks.
The whole list is fitted once, after it has loaded. See
[Peak Lists and Output Files](peaklistformats.md).

You can also load a list as the window opens, by passing it to any 2D routine:

```julia
relaxation2d("11/pdata/1"; peaklist="out/peaklist.csv")
```

On Linux, the **Load peak list** button doesn't open a file dialog, because the GTK file
dialog can crash Julia there; it prints a warning instead. Pass `peaklist` when you start
the analysis. The crash comes from GTK settings schemas older than the GTK that Julia
bundles, so updating your system's GTK 3 package (`libgtk-3-common` on Debian or Ubuntu,
`gtk3` on Fedora) may fix the dialog.

## Summary plots

[`summaryplot`](summary.md) plots a fitted parameter against residue number,
from a live experiment or one or more saved `results.csv` files. See the
[Summary Plots](summary.md) page for full details and examples.

## Adjusting the Fitting Region

The X and Y radius sliders in the peak info panel set the default size of the region
around each peak used for lineshape fitting. Smaller radii are appropriate for crowded
spectra; larger radii improve the fit for broad peaks.

To give one peak radii of its own, select it and press `Shift` with the arrow keys:
left and right narrow and widen the x radius, down and up the y radius. A peak with its own
radii keeps them when you change the sliders, and the info panel shows them. Press `=` to
return it to the defaults.
