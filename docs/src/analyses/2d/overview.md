# 2D Experiment Analysis

NMRAnalysis.jl provides interactive GUIs for analysing 2D NMR experiments, including
relaxation, exchange, and NOE measurements. All functions follow the same pattern: they
load one or more processed 2D spectra, open an interactive window for peak picking and
fitting, and export results to a folder of your choice.

The [experiment-specific pages](fit.md) describe the available functions
and the theory behind each analysis. This page covers the shared GUI features common to
all of them.

<!-- screenshot: assets/2d-gui-overview.png -->

## Adding and Managing Peaks

Peaks are picked interactively using the mouse and keyboard. Move the cursor over a peak
in the contour plot to work with it.

| Action | Key / Button |
|--------|--------------|
| Add peak at cursor position | `A` |
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
| `results.csv` | One row per peak: positions (δ₁, δ₂), linewidths (R2x, R2y), per-plane amplitudes, and any derived experiment parameters (relaxation rates, NOE values, …), each with uncertainties |
| `summary.pdf` | Summary plot of the primary fitted parameter against residue number (or atom for methyl/non-backbone experiments) |
| `peak_LABEL.pdf` | Per-peak publication-quality fit plot for each labelled peak |
| `cluster_LABEL.pdf` | Zoomed 2D contour plot (first plane) with fitted lineshapes for each group of overlapping peaks |

`results.csv` has experiment metadata in `#`-comment lines above an ordinary
header row, so it opens directly in spreadsheets and `pandas`. Existing files are
backed up with an `.old` extension before being overwritten. See
[Peak Lists and Output Files](peaklistformats.md) for the full column description.

## Loading and Resuming Analysis

The **Load peak list** button restores peak positions and labels from a saved
`results.csv` (or a simple `label x y` text file), so you can resume work later or
seed a new analysis from existing positions. Only the `label`, `x` and `y` columns
are read — see [Peak Lists and Output Files](peaklistformats.md).

## Summary plots

[`summaryplot`](summary.md) plots a fitted parameter against residue number,
from a live experiment or one or more saved `results.csv` files. See the
[Summary Plots](summary.md) page for full details and examples.

## Adjusting the Fitting Region

The X and Y radius sliders in the peak info panel control the size of the region around
each peak used for lineshape fitting. Smaller radii are appropriate for crowded spectra;
larger radii improve the fit for broad peaks.
