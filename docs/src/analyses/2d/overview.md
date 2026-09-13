# 2D Experiment Analysis

NMRAnalysis.jl provides interactive GUIs for analysing 2D NMR experiments, including
relaxation, exchange, and NOE measurements. All functions follow the same pattern: they
load one or more processed 2D spectra, open an interactive window for peak picking and
fitting, and export results to a folder of your choice.

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
| Track peak from cursor position | `T` |
| Delete the selected peak | `D` or **Delete peak** button |
| Rename the selected peak | `R` or **Rename peak** button |
| Navigate to previous spectrum slice | `←` or **←** button |
| Navigate to next spectrum slice | `→` or **→** button |
| Raise contour base level | `↑` or **contour ↑** button |
| Lower contour base level | `↓` or **contour ↓** button |
| Reset axis zoom | **reset zoom** button |
| Show or hide the fitted lineshape overlay | **Fitting** toggle |
| Open a summary plot of the current results | **Summary plot** button (enabled once peaks are present) |
| Load a previously saved peak list | **Load peak list** button |
| Save all results to a folder | **Save to folder** button |
| Close the GUI window | **Quit** button |

Peak lineshapes are fitted in real time as you add or move peaks. The right panel shows
cross-sections (or a model fit plot, for relaxation-type experiments) for the currently
selected peak.

## Visual Feedback

The window background changes colour to indicate the current interaction mode:

| Background | Mode |
|------------|------|
| White | Normal |
| Salmon / orange | Fitting in progress (save operation) |
| Light blue | Renaming a peak |
| Pale green | Moving a peak |

Peak markers are colour-coded:

| Colour | Meaning |
|--------|---------|
| Blue | Unmodified peak |
| Red | Manually moved or adjusted peak |
| Green | Currently selected peak |

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
by plane, the whole trajectory is restored. See
[Peak Lists and Output Files](peaklistformats.md).

## Summary plots

[`summaryplot`](summary.md) plots a fitted parameter against residue number,
from a live experiment or one or more saved `results.csv` files. See the
[Summary Plots](summary.md) page for full details and examples.

## Adjusting the Fitting Region

The X and Y radius sliders in the peak info panel control the size of the region around
each peak used for lineshape fitting. Smaller radii are appropriate for crowded spectra;
larger radii improve the fit for broad peaks.
