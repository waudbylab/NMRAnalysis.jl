# Quick Start

## Installation

1. **Install Julia** from [julialang.org](https://julialang.org/install/) if you haven't already.

2. **Add NMRAnalysis.jl** from the Julia package manager:
   ```
   julia> ]
   pkg> add NMRAnalysis
   ```
   Julia will download and compile NMRAnalysis and its dependencies. This takes a few
   minutes the first time.

3. **Load the package**:
   ```julia
   using NMRAnalysis
   ```
   A startup message will list all available analysis functions.

## Basic Workflow

Most analyses follow the same pattern: navigate to your data directory, then call the
appropriate function.

```julia
using NMRAnalysis

cd("/path/to/your/nmr/data")
```

Get help on any function with `?`:
```julia
?relaxation1d
?relaxation2d
?r1rho
```

## Examples

### 1D Relaxation (R₁ / R₂)

```julia
relaxation1d("11")                          # delays from the vdlist or annotations
relaxation1d("11"; tau=[0.01, 0.05, 0.1])   # or given explicitly
```

You are asked about anything that cannot be read from the experiment, and a window then
opens in which you position the integration region and watch the fit. The same pattern
applies to `diffusion1d`, `tract`, `calibration1d` and `kinetics1d`. See
[1D Experiments](analyses/1d/overview.md).

### ¹⁹F R1ρ Relaxation Dispersion

```julia
r1rho()                       # select experiment folder interactively
r1rho("11")                   # specify experiment directly
r1rho("11"; minvSL=500)       # filter low spin-lock powers
```

See the [R1ρ Tutorial](tutorials/r1rho.md) for a step-by-step guide.

### 2D Relaxation (T1 / T2)

```julia
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1"],
    [0.010, 0.030, 0.060, 0.100]   # delays in seconds
)
```

An interactive graphical window opens for peak picking and fitting. See
[Relaxation (T1/T2)](analyses/2d/relaxation.md) for details.

## Next Steps

- **[1D Experiments](analyses/1d/overview.md)** - the shared analysis window, and a guide for each 1D experiment type
- **[2D Experiments](analyses/2d/overview.md)** — interactive GUI reference and per-experiment guides
- **[Tutorials](tutorials/r1rho.md)** — step-by-step worked examples

## Getting Help

- Use `?function_name` in the Julia REPL for built-in help
- Report issues and suggest features at [github.com/waudbylab/NMRAnalysis.jl](https://github.com/waudbylab/NMRAnalysis.jl/issues)
