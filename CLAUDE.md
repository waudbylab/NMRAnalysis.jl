# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NMRAnalysis.jl is a Julia package for analysis of NMR experiments, specifically diffusion and relaxation experiments. The package provides both 1D and 2D experiment analysis capabilities with interactive GUI interfaces.

## Development Commands

### Package Management
- `julia --project=.` - Start Julia with the project environment
- `julia --project=. -e "using Pkg; Pkg.instantiate()"` - Install dependencies
- `julia --project=. -e "using Pkg; Pkg.test()"` - Run tests (currently not working)

### Documentation
- `julia --project=docs docs/make.jl` - Build documentation locally
- Documentation is built with Documenter.jl and deployed to GitHub Pages

### Code Formatting
- Uses JuliaFormatter with "yas" style (configured in `.JuliaFormatter.toml`)
- `julia --project=. -e "using JuliaFormatter; format(\".\")"` - Format all code

## Architecture

### Main Module Structure
- **Main module**: `src/NMRAnalysis.jl` - Entry point that includes and re-exports functionality
- **GUI2D module**: `src/gui2d/` - 2D experiment analysis with interactive plotting
- **R1rho module**: `src/R1rho/` - R1ρ relaxation dispersion experiments
- **1D analysis functions**: Individual files for diffusion, TRACT, viscosity analysis

### Core Components

#### 1D Experiment Analysis
- `diffusion()` - Diffusion coefficient analysis
- `tract()` - TRACT (Temperature-Ramped Analysis of Conformational Transitions)
- `r1rho()` - R1ρ relaxation dispersion analysis
- `viscosity` - Viscosity calculations

#### 2D Experiment Analysis
All 2D functions are provided by the GUI2D module:
- `intensities2d()` - Basic intensity measurements
- `relaxation2d()` - Relaxation parameter fitting
- `recovery2d()` - Recovery experiments
- `modelfit2d()` - Model fitting interface
- `hetnoe2d()` - Heteronuclear NOE experiments
- `cest2d()` - CEST (Chemical Exchange Saturation Transfer)
- `cpmg2d()` - CPMG (Carr-Purcell-Meiboom-Gill) experiments
- `pre2d()` - PRE (Paramagnetic Relaxation Enhancement)

### Key Dependencies
- **NMRTools.jl**: Core NMR data handling and processing
- **Makie.jl** (CairoMakie/GLMakie): Interactive plotting and GUI interfaces
- **LsqFit.jl**: Nonlinear least squares fitting
- **Measurements.jl/MonteCarloMeasurements.jl**: Uncertainty propagation
- **NativeFileDialog.jl**: File selection dialogs

### GUI Architecture
- Both R1rho and GUI2D modules use Makie for interactive GUIs
- State management pattern with modular event handling
- Mouse and keyboard interaction systems
- Real-time parameter fitting and visualization

### Module Organization
- `GUI2D/` contains ~25 files organizing different aspects of 2D analysis
- `R1rho/` contains ~10 files for R1ρ-specific analysis
- Each module has its own state management, GUI, fitting, and experiment handling

### Data Flow
1. File selection via `select_expts()` or file dialogs
2. Data processing and experiment setup
3. Interactive GUI for peak picking, parameter adjustment
4. Real-time fitting and result visualization
5. Export capabilities for results

## Common Patterns

### Function Entry Points
Most analysis functions follow the pattern:
```julia
function_name(directory_path=""; options...)
function_name(filenames::Vector{String}; options...)
```

### GUI Initialization
1. Activate appropriate Makie backend (GLMakie for interactive, CairoMakie for static)
2. Process experiment files into dataset
3. Initialize state object
4. Launch GUI with `gui!(state)`

### Error Handling
- Extensive validation of input parameters
- Graceful handling of file selection cancellation
- Informative error messages for common issues

### Naming Conventions
- Follow the Julia style guide: function names should be lowercase without underscores (e.g. `plot(results)`, not `plot_result(results)`)
- Do NOT prefix private/internal functions with `_`. Just use unexported names.
- Avoid underscores in function names entirely — use concatenated lowercase words (e.g. `combineplots`, `defaultparams`, `formatvalue`)