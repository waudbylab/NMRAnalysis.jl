# NMRAnalysis.jl

NMRAnalysis.jl is a Julia package for interactive analysis of biomolecular NMR experiments,
with a particular focus on relaxation, diffusion, exchange, and protein dynamics measurements.
It provides both 1D command-line workflows and 2D interactive graphical interfaces, and is
designed to make common NMR analyses straightforward without sacrificing flexibility.

!!! note "Active development"
    NMRAnalysis.jl is under active development. Features and API may change as the
    package evolves.

## What's available

### 1D Analyses
Command-line tools for routine 1D experiments:

| Function | Experiment |
|----------|-----------|
| [`diffusion()`](analyses/diffusion.md) | DOSY / diffusion coefficient measurement |
| [`relaxation()`](analyses/relaxation.md) | R1 and R2 relaxation (exponential or inversion-recovery fit) |
| [`tract()`](analyses/tract.md) | TRACT experiment for rotational correlation time |
| [`r1rho()`](analyses/r1rho.md) | ¹⁹F R1ρ relaxation dispersion |

### 2D Analyses
Interactive graphical interfaces for 2D and pseudo-3D experiments:

| Function | Experiment |
|----------|-----------|
| [`fit2d()`](analyses/2d/fit.md) | Peak fitting: positions, linewidths, and amplitudes |
| [`peaktrack2d()`](analyses/2d/peaktracking.md) | Tracking peak positions across a series of 2D spectra |
| [`titration2d()`](analyses/2d/titration.md) | Binding isotherms from a titration series |
| [`relaxation2d()`](analyses/2d/relaxation.md) | R1 / R2 relaxation from a series of 2D spectra |
| [`recovery2d()`](analyses/2d/magnetisationrecovery.md) | Inversion or saturation recovery |
| [`modelfit2d()`](analyses/2d/modelfit.md) | Custom model fitting |
| [`hetnoe2d()`](analyses/2d/hetnoe.md) | Heteronuclear NOE |
| [`cpmg2d()`](analyses/2d/cpmg.md) | CPMG relaxation dispersion |
| [`cest2d()`](analyses/2d/cest.md) | CEST (Chemical Exchange Saturation Transfer) |
| [`pre2d()`](analyses/2d/pre.md) | Paramagnetic relaxation enhancement |
| [`ccr2d()`](analyses/2d/ccr.md) | Cross-correlated relaxation |
| [`methylccr2d()`](analyses/2d/methylccr.md) | Methyl CCR: S²τc from buildup/decay series |
| [`rdc2d()`](analyses/2d/rdc.md) | Residual dipolar couplings |

All 2D functions share the same interactive GUI — see the [2D Overview](analyses/2d/overview.md) for
a guide to the interface and common workflow.

## Ecosystem

NMRAnalysis.jl is part of a suite of Julia packages for NMR data handling developed
by the [Waudby lab](https://waudbylab.org):

- **[NMRTools.jl](https://waudbylab.org/NMRTools.jl)** — the foundation for
  NMR data import and processing in Julia. NMRAnalysis.jl is built on top of NMRTools
  for all data loading, axis handling, and spectral processing.

- **[NMRScreen.jl](https://waudbylab.org/NMRScreen.jl)** — tools for
  fragment and ligand screening by NMR, including automated analysis of large
  compound libraries.

## Contributing

NMRAnalysis.jl is developed and maintained by the
[Waudby lab](https://waudbylab.org) at University College London.
Contributions are warmly welcomed — whether that's bug reports, new analysis
routines, documentation improvements, or example datasets. Please open an issue
or pull request on [GitHub](https://github.com/waudbylab/NMRAnalysis.jl).
