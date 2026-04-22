# NMRAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://waudbylab.github.io/NMRAnalysis.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://waudbylab.github.io/NMRAnalysis.jl/dev)
[![CI](https://github.com/waudbylab/NMRAnalysis.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/waudbylab/NMRAnalysis.jl/actions/workflows/Runtests.yml)
[![DOI](https://zenodo.org/badge/665599660.svg)](https://doi.org/10.5281/zenodo.15667046)

![NMRAnalysis logo](logo.png)

NMRAnalysis.jl is a library for analysis of NMR experiments. This includes 1D experiments:
- diffusion
- relaxation
- TRACT
- R1rho relaxation dispersion

and 2D experiments:
- lineshape fitting (for measurement of intensities and linewidths)
- relaxation
- heteronuclear NOEs
- CPMG relaxation dispersion
- CEST


## System requirements

- Julia 1.10 (Julia 1.12 recommended)
- Developed and tested using OS X, additional testing carried out on linux systems


## Installation

1. If you don't already have Julia installed, **download Julia** from [https://julialang.org/install/](https://julialang.org/install/) and install according to the instructions on this page.

2. **Verify installation** by opening a terminal/command prompt and typing:
   ```
   julia
   ```
   You should see the Julia REPL (interactive prompt) with version information.

3. Once Julia is installed, you need to **add the NMRAnalysis package**:

   1. Enter package mode by pressing `]` (you'll see the prompt change to `pkg>`)
   2. Add the package by typing:
      ```
      add NMRAnalysis
      ```
   3. Wait for installation - Julia will automatically download and install NMRAnalysis and all its dependencies
   4. Exit package mode by pressing Backspace

4. **Activate the installation** by typing:
   ```julia
   using NMRAnalysis
   ```

   You should see an information message listing the available analysis routines.

> **NOTE**: The first time you install NMRAnalysis, it may take several minutes to download and compile all dependencies. This is normal and only happens once.


## Demo

Tutorials are available in the online documentation, including example data:
 - 19F R1rho relaxation dispersion analysis of ligand binding: [https://waudbylab.org/NMRAnalysis.jl/stable/tutorials/r1rho/](https://waudbylab.org/NMRAnalysis.jl/stable/tutorials/r1rho/)

## Instructions for use

1. **Navigate to your data directory**:
   ```julia
   cd("/path/to/your/nmr/data")  # Replace with your actual data path
   ```
   For example: `cd("/Users/chris/NMR/my_experiment")`

2. **Load the package**:
   ```julia
   using NMRAnalysis
   ```

3. **Get help** on any function:
   ```julia
   ?diffusion    # Shows help for diffusion analysis
   ?r1rho        # Shows help for R1ρ analysis
   ```

Some examples of analysing **1D experiments**:

```julia
using NMRAnalysis

# Diffusion analysis - analyzes DOSY experiments
diffusion("106")  # Analyze experiment in folder "106"

# TRACT analysis - for rotational correlation times
tract()           # Prompt to select experiment folders

# R1ρ relaxation dispersion
r1rho()                                   # Show file selection dialog
r1rho(["11", "12"])                       # Analyze experiments 11 and 12
r1rho(["11", "12"], minvSL=500)           # Filter low spin-lock strengths
```

**2D experiments** can be analysed with the interactive graphical interface:

```julia
# Relaxation analysis (T1, T2)
relaxation2d(
    "expno",            # Processed spectra as pseudo-3d
    [0.01, 0.03, 0.05]  # Relaxation delays (s)
)

# Heteronuclear NOE analysis
hetnoe2d(
    ["reference/pdata/1", "saturated/pdata/1"],  # Reference and saturated spectra
    [false, true]                                # Saturation states
)
```

### File formats

NMRAnalysis works with standard Bruker data formats:

```julia
# Single experiment (TopSpin experiment number)
diffusion("106")

# Multiple experiments
r1rho(["11", "12", "13"])

# Processed data directories
hetnoe2d(["reference/pdata/1", "saturated/pdata/1"], [false, true])

# Full paths (if data is elsewhere)
diffusion("/Users/chris/NMR/project_data/106")
```

## Documentation

Access the online documentation at [https://waudbylab.org/NMRAnalysis.jl/stable/](https://waudbylab.org/NMRAnalysis.jl/stable/)

