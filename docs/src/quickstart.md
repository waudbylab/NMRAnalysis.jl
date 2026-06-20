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
?relaxation2d
?r1rho
```

## 1D Experiments

### Diffusion (DOSY)

```julia
diffusion("106")              # analyse experiment 106
diffusion()                   # select folder interactively
```

See [Diffusion](@ref "1D Diffusion Analysis") for a full walkthrough.

### R1 / R2 Relaxation

```julia
relaxation("5")               # exponential decay or inversion recovery
```

See [Relaxation](@ref "1D Relaxation Analysis") for details.

### TRACT

```julia
tract()                       # select TROSY and anti-TROSY experiments interactively
```

### ¹⁹F R1ρ Relaxation Dispersion

```julia
r1rho()                       # select experiment folder interactively
r1rho("11")                   # specify experiment directly
r1rho("11"; minvSL=500)       # filter low spin-lock powers
```

See the [R1ρ Tutorial](@ref "19F R1ρ acquisition & analysis") for a step-by-step guide.

## 2D Experiments

All 2D analysis functions open an interactive graphical window. Use the mouse to pick
peaks and the keyboard shortcuts described in the [2D Overview](@ref "Overview") to
navigate and manage peaks.

### Peak fitting (positions, linewidths, amplitudes)

```julia
fit2d("expno/pdata/1")                              # single spectrum
fit2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"])  # series of spectra
```

### Relaxation (T1 / T2)

```julia
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1"],
    [0.010, 0.030, 0.060, 0.100]   # delays in seconds
)
```

### Heteronuclear NOE

```julia
hetnoe2d(
    ["ref/pdata/1", "sat/pdata/1"],
    [false, true]
)
```

### CPMG Relaxation Dispersion

```julia
cpmg2d("11/pdata/1"; Trelax=0.04, vCPMG=[0, 25, 50, 100, 200, 500])
```

### CEST

```julia
cest2d("11/pdata/1"; B1=15, Tsat=0.3)
```

## Next Steps

- **[1D Analysis](@ref "Diffusion")** — detailed guides for each 1D experiment type
- **[2D Analysis](@ref "Overview")** — interactive GUI reference and per-experiment guides
- **[Tutorials](@ref "19F R1ρ acquisition & analysis")** — step-by-step worked examples

## Getting Help

- Use `?function_name` in the Julia REPL for built-in help
- Report issues and suggest features at [github.com/waudbylab/NMRAnalysis.jl](https://github.com/waudbylab/NMRAnalysis.jl/issues)
