# NMRAnalysis.jl

NMRAnalysis.jl is a Julia package for interactive analysis of biomolecular NMR experiments,
with a particular focus on relaxation, diffusion, exchange, and protein dynamics measurements.
It provides both 1D command-line workflows and 2D interactive graphical interfaces, and is
designed to make common NMR analyses straightforward without sacrificing flexibility.

!!! note "Active development"
    NMRAnalysis.jl is actively developed and extended. The core analysis functions are
    stable and in routine use, but the API may evolve as new features are added.

## What's available

### 1D Analyses
Command-line tools for routine 1D experiments:

| Function | Experiment |
|----------|-----------|
| [`diffusion()`](@ref) | DOSY / diffusion coefficient measurement |
| [`relaxation()`](@ref) | R1 and R2 relaxation (exponential or inversion-recovery fit) |
| [`tract()`](@ref) | TRACT experiment for rotational correlation time |
| [`r1rho()`](@ref) | ¹⁹F R1ρ relaxation dispersion |

### 2D Analyses
Interactive graphical interfaces for 2D and pseudo-3D experiments:

| Function | Experiment |
|----------|-----------|
| [`fit2d()`](@ref) | Peak fitting: positions, linewidths, and amplitudes |
| [`relaxation2d()`](@ref) | R1 / R2 relaxation from a series of 2D spectra |
| [`recovery2d()`](@ref) | Inversion or saturation recovery |
| [`modelfit2d()`](@ref) | Custom model fitting |
| [`hetnoe2d()`](@ref) | Heteronuclear NOE |
| [`cpmg2d()`](@ref) | CPMG relaxation dispersion |
| [`cest2d()`](@ref) | CEST (Chemical Exchange Saturation Transfer) |
| [`pre2d()`](@ref) | Paramagnetic relaxation enhancement |
| [`ccr2d()`](@ref) | Cross-correlated relaxation |

All 2D functions share the same interactive GUI — see [Overview](@ref "Overview") for
a guide to the interface and common workflow.

## Ecosystem

NMRAnalysis.jl is part of a suite of Julia packages for NMR data handling developed
by the [Waudby lab](https://waudbylab.org):

- **[NMRTools.jl](https://github.com/waudbylab/NMRTools.jl)** — the foundation for
  NMR data import and processing in Julia. NMRAnalysis.jl is built on top of NMRTools
  for all data loading, axis handling, and spectral processing.

- **[NMRScreen.jl](https://github.com/waudbylab/NMRScreen.jl)** — tools for
  fragment and ligand screening by NMR, including automated analysis of large
  compound libraries.

## Contributing

NMRAnalysis.jl is developed and maintained by the
[Waudby lab](https://waudbylab.org) at University College London.
Contributions are warmly welcomed — whether that's bug reports, new analysis
routines, documentation improvements, or example datasets. Please open an issue
or pull request on [GitHub](https://github.com/waudbylab/NMRAnalysis.jl).
