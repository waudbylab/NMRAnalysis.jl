# NMRAnalysis.jl

Fit relaxation, diffusion, exchange and titration experiments from Bruker data. Each
experiment has one function: call it with your data, pick peaks or regions in the window
that opens, and save the results.

```julia
using NMRAnalysis
relaxation2d(example("relaxation2d_asyn"); cycletime=0.02809)
```

![The relaxation2d window: an ¹H,¹⁵N spectrum with picked peaks, and the fitted decay of the selected peak](relaxation2d.png)

To install the package and run your first analysis, see the [Quick start](quickstart.md).

!!! note "Active development"
    NMRAnalysis.jl is under active development, and its interface may change.

## Contributing

NMRAnalysis.jl is developed by the [Waudby lab](https://waudbylab.org) at University College
London. Report bugs, request features or contribute example datasets on
[GitHub](https://github.com/waudbylab/NMRAnalysis.jl).
