# 1D Experiments

The 1D routines analyse a series of 1D spectra recorded as a single pseudo-2D experiment:
a relaxation decay, a diffusion ramp, a nutation calibration, a reaction followed over
time. Each is a single function that you call from the Julia prompt with the experiment
folder. Every one of them works the same way:

1. It reads the experiment and asks about anything it needs and cannot work out for
   itself.
2. It opens a window showing your spectra, with an integration region you position over
   the signal you care about.
3. It fits as you work, so you can see immediately whether the region you have chosen
   gives a sensible result.
4. It saves the fit, the numbers and a results file to a folder when you press **Save**.

<!-- TODO screenshot (1): docs/src/assets/analysis1d-window.png
     The relaxation1d window in a typical state: spectra overlaid on the left with one
     shaded integration region and the noise marker, the fitted decay below it, and the
     results panel on the right. This is the page's lead image. -->

## Which routine to use

| Function | Experiment | What you get |
|---|---|---|
| [`relaxation1d`](relaxation.md) | T₁ or T₂ decay, inversion recovery | relaxation rate `R` |
| [`diffusion1d`](diffusion.md) | diffusion (DOSY) gradient ramp | diffusion coefficient `D`, hydrodynamic radius `rH` |
| [`tract`](tract.md) | TROSY and anti-TROSY pair | rotational correlation time `τc` |
| [`calibration1d`](calibration.md) | nutation (pulse-length) calibration | 90° pulse length, ν₁, B₁ inhomogeneity |
| [`kinetics1d`](kinetics.md) | any time series | intensity of each region against time |

Chemical exchange experiments have their own routines: see
[R1ρ Relaxation Dispersion](../r1rho.md) and
[Exchange Fitting](../exchange1d/overview.md).

## Your first analysis

Start Julia, move to your data, and call the routine for your experiment:

```julia
using NMRAnalysis

cd("/Users/chris/NMR/my_dataset")
results = relaxation1d("11")
```

The routine reads experiment 11, asks about anything it needs (see
[Where the parameters come from](#Where-the-parameters-come-from) below), and opens the
analysis window. Position the region over your signal, check the fit, press **Save**, and
close the window. `results` then holds the fitted parameters.

## Where the parameters come from

Every analysis needs to know something beyond the spectra themselves: the relaxation
delays, the gradient strengths, which model to fit. Each routine looks for each of these
in the same order and stops at the first one that answers:

| Order | Source | Example |
|---|---|---|
| 1 | An argument you passed | `relaxation1d("11"; tau=[0.01, 0.05, 0.1])` |
| 2 | An annotation in the pulse sequence | `relaxation.duration` |
| 3 | A Bruker acquisition parameter | the `vdlist`, `p30`, `d20` |
| 4 | A question, asked before the window opens | *"Enter 8 relaxation delays..."* |

Nothing is required to be annotated. Annotations, where a pulse sequence carries them,
simply save you from being asked. Anything you pass as an argument is taken as given and
is not asked about, so a working analysis can always be turned into a script that repeats
itself without any typing.

!!! note "Confirming what was found"
    [`diffusion1d`](diffusion.md) is the exception that also offers what it *did* find for
    confirmation, because its parameters are easy to get wrong and expensive to get wrong
    silently. Press enter to accept each one, or type a replacement.

In a script, or anywhere else that is not an interactive Julia session, no question is
ever asked. A value that has a sensible default takes it; one that does not raises an
error naming the argument to pass instead. Force either behaviour with the `prompt`
keyword:

```julia
relaxation1d("11"; prompt=false)   # never ask; error if something is missing
```

## Using the analysis window

The window has three parts. On the left, your spectra are overlaid in grey with the
current one picked out in colour, and below them the fit for whichever region is selected.
On the right are the controls and a panel giving the fitted numbers.

The **shaded regions** are what get integrated. The selected one is orange and the rest
are blue; whichever the cursor is over becomes the selected one. The **purple marker** is
where the noise is measured. It always uses the same width as the region being integrated,
which is what makes it a direct estimate of that region's uncertainty, so put it somewhere
with no signal in any spectrum.

### Mouse

| To do this | Do that |
|---|---|
| Select a region | Move the cursor over it |
| Move a region | Drag its middle |
| Resize a region | Drag either edge, or hold shift and scroll |
| Move the noise marker | Drag it |
| Zoom in on part of the spectrum | Drag a rectangle, or scroll |

### Keyboard

| Key | Action |
|---|---|
| `A` | Add a region. Tap for a default width, or press, drag and release to size it |
| `R` | Rename the selected region |
| `D` | Delete the selected region |
| `←` `→` | Step through the spectra |
| `↑` `↓` | Double or halve the vertical scale |

The background turns yellow while you are dragging out a new region, and blue while you
are typing a name. Press escape to cancel a rename.

### Buttons

| Button | Action |
|---|---|
| **×2**, **÷2**, **reset zoom** | Vertical scale, and back to the full spectrum |
| **←**, **→**, slider | Step through the spectra |
| **(R)ename**, **(D)elete** | The same as the keys, for the selected region |
| **Fitting** | Turn fitting off to look at the raw integrals alone |
| **Save** | Write every output file to the output folder |
| **Load** | Restore the regions saved in that folder by an earlier session |

The text box beside **Save** names the output folder, relative to your working directory.
It defaults to `out`.

## Suggested workflow

1. Step to a spectrum where the signal is clear. For an inversion recovery, that is
   usually the last one.
2. Drag the region over the signal you want, and resize it to cover the whole peak.
3. Drag the noise marker to an empty part of the spectrum.
4. Check the fit below. The points should sit on the line, and the error bars should look
   like the scatter.
5. Add further regions with `A` if you want more than one signal, and name them with `R`.
6. Press **Save**, then close the window.

!!! tip "Choosing a region"
    Integrate the whole peak, including its wings, rather than just the top of it. For
    TRACT, integrate the entire amide envelope, which is what the default 7.5 to 9.5 ppm
    region does.

## What gets saved

**Save** writes four kinds of file into the output folder, backing up an existing
`results.csv` as `results.csv.bak` first:

| File | Contents |
|---|---|
| `results.csv` | One row per region (and per series, where an experiment has more than one), with every fitted and derived parameter and its uncertainty |
| `summary.txt` | The same numbers laid out to be read, with the region and noise positions |
| `fit.pdf` | The fit for every region on one set of axes |
| `fit_<region>.pdf` | The fit for each region on its own |

`results.csv` carries the experiment details as `#` comment lines above an ordinary header
row, so it opens directly in a spreadsheet or `pandas`. It is also what **Load** reads: it
restores the regions and the noise position, so a session with several named regions can
be picked up again later rather than picked out by hand.

## Repeating an analysis

Each routine returns a vector of results, one per region and series, whether you used the
window or not:

```julia
results = relaxation1d("11")
param(results[1], :R)     # the fitted rate, with its uncertainty
```

`param` reads both the fitted parameters and anything derived from them, so
`param(r, :τc)` and `param(r, :rH)` work the same way.

To repeat an analysis without the window at all, pass the region you settled on as an
`integration` triple. This is the same triple the exchange routines use, so a region
chosen once can serve a whole set of experiments:

```julia
relaxation1d("11"; integration=(peakppm=8.2, noiseppm=-1.0, ppmwidth=0.5))
```

[`pickregion`](@ref) opens the region-picking window on its own, with no fitting, and
returns that triple.

## Uncertainties

Intensities are divided by the spectrum's own noise level when the data is loaded, so the
numbers in the window are in units of signal to noise. The uncertainty on an *integrated*
region is measured separately, and directly: the same width of spectrum is integrated at
the noise marker in every spectrum of the series, and the standard deviation of those
integrals is taken. That is why the noise marker matters, and why it tracks the width of
the region it is estimating noise for.

This assumes a flat baseline under the noise marker. If the fitted uncertainties look
implausible, that is the first thing to check.

## Automatic dispatch

Where a pulse sequence carries the annotations for it, [`analyse`](../analyse.md) will
recognise the experiment and call the right routine for you:

```julia
analyse("11")
```

This works today for nutation calibrations and for R1 and R2 relaxation experiments.
Diffusion and TRACT are called directly.
