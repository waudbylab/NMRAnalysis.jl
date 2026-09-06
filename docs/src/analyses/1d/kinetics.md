# Kinetics

`kinetics1d` follows the intensity of one or more regions of a spectrum against time. It
is the general-purpose routine for a reaction, an exchange process, or anything else
recorded as a series of 1D spectra: rather than fitting a particular model, it gives you
the intensity trace for each region you mark, which you can then fit or plot yourself.

## Running the analysis

`kinetics1d` needs the times, because there is no acquisition parameter that reliably
holds them:

```julia
using NMRAnalysis

times = 0:120:7200                      # seconds, one per spectrum
results = kinetics1d("21", times)
```

The [analysis window](overview.md#Using-the-analysis-window) opens with one region on the
tallest peak. Add a region for each signal you want to follow by pressing `A`, name it
with `R`, and press **Save**.

<!-- TODO screenshot (6): docs/src/assets/kinetics1d-window.png
     kinetics1d with three named regions over separate signals, their intensity traces
     plotted below in their own colours. There is no existing screenshot to replace. -->

## Several runs

Where a series contains repeats, tag each spectrum with `run` and each run is treated as
its own series and plotted in its own colour:

```julia
kinetics1d("21", times; run=[1, 1, 1, 2, 2, 2])
```

## Fitting a model

By default nothing is fitted: the deliverable is the intensity against time. Pass a
`SeriesModel` to fit one, for instance a single exponential decay:

```julia
kinetics1d("21", times; model=ExponentialModel())
```

## Reading the traces

The returned results carry the trace itself, one entry per region and run:

```julia
results = kinetics1d("21", times)
r = results[1]
r.x        # times
r.y        # integrated intensities, with uncertainties
```

!!! note "The traces are not written to `results.csv` yet"
    `results.csv` holds one row per region and run with the fitted parameters, which for
    an unfitted kinetics series means the region bounds alone. Take the traces from the
    returned results, as above, until the output format for a series-valued result is
    settled.
