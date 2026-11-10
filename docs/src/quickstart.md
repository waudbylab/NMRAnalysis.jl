# Quick start

## Install

Install [Julia](https://julialang.org/install/), then add the package:

```
julia> ]
pkg> add NMRAnalysis
```

The first install and load take a few minutes while Julia compiles the package.

## Analyse a 2D relaxation series

This example uses a ¹⁵N R₂ series recorded on α-synuclein. `example` downloads it the first
time you use it.

```julia
using NMRAnalysis
expt = relaxation2d(example("relaxation2d_asyn"); cycletime=0.02809)
```

The relaxation delays come from the experiment's `vclist`, multiplied by `cycletime`, the
duration of one loop in seconds. Without `cycletime`, you're asked for it.

In the window that opens:

1. Move the cursor over a peak and press `A`. The peak is fitted, and its decay appears on
   the right.
2. Press `R` to rename the peak, for example `V40`, so that plots use the residue number.
3. Add the other peaks you want.
4. Press **Save**. The results go into the folder named in the box beside it, `out` by
   default.
5. Press `Q` to close the window.

```@setup quickstart
using NMRAnalysis
expt = NMRAnalysis.screenshot("relaxation2d.png"; peak="V40") do
    reopen(exampleanalysis("relaxation2d_asyn"))
end
```

![The relaxation2d window with the example's peaks picked](relaxation2d.png)

When the window closes, `relaxation2d` returns the analysis. `results(expt)` gives one row
per peak:

```@example quickstart
Main.resultstable(results(expt)) # hide
```

To open the analysis shown above, with its peaks already picked, run:

```julia
reopen(exampleanalysis("relaxation2d_asyn"))
```

`reopen("out")` reopens your own analysis in the same way.

## Analyse a 1D experiment

1D experiments work in the same way, with integration regions in place of peaks. This
example measures the rotational correlation time of calreticulin from a TRACT pair:

```julia
trosy, antitrosy = example("tract_calr")
tract(trosy, antitrosy)
```

```@setup quickstart
NMRAnalysis.screenshot("tract.png") do
    tract(example("tract_calr")...)
end
```

![The TRACT window, with the TROSY and anti-TROSY decays fitted](tract.png)

## Use your own data

Change to the folder holding your experiments, and give experiment numbers:

```julia
cd("/path/to/data")
relaxation2d(11)
```

Each function reads what it can from the acquisition parameters and the pulse-sequence
annotations, and asks for anything else. Type `?relaxation2d` for the options.
