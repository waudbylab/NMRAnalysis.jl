# Kinetics

`kinetics1d` follows the integral of one or more regions over time, for a reaction or any
other process recorded as a series of 1D spectra. It fits nothing unless you give a model:
the result is the intensity of each region against time.

## Usage

Give the time of each spectrum, in seconds:

```julia
using NMRAnalysis
res = kinetics1d("21"; T=0:120:7200)
```

Add a region for each signal with `A`, and name it with `R`.

| Keyword | Meaning |
|---|---|
| `T` | Time of each spectrum, s |
| `run` | Which run each spectrum belongs to, where the series holds repeats, e.g. `[1, 1, 1, 2, 2, 2]` |
| `model` | A model to fit, such as `ExponentialModel()` |

## Output

`series.csv` holds the integral of each region at each time:

```
source,label,time (s),I,I_err,I_fit
21/pdata/1,reactant,0.0,9421.3,12.4,NA
21/pdata/1,reactant,60.0,7724.5,12.4,NA
21/pdata/1,product,0.0,0.0,12.4,NA
21/pdata/1,product,60.0,1702.1,12.4,NA
```

The returned results hold the same: `res[1].x` is the times, and `res[1].y` the integrals,
with uncertainties.
