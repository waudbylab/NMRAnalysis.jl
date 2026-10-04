"""
    NMRAnalysis

Interactive analysis of NMR relaxation, diffusion, exchange and 2D peak series. Set your
working directory to the folder holding your experiments, call one of the routines below,
and use `?name` for help on any of them. Bruker experiment numbers work wherever a path
does, and any parameter that is not given is read from the pulse-sequence annotations or
acquisition parameters where possible, and otherwise asked for.

# Automatic analysis

- `analyse(filename)`, `analyse([filenames...])`: choose a routine from the annotations

# 1D experiments

Each opens a window to pick integration and noise regions, fits live, and returns the
results when the window closes.

- `relaxation1d(filename)`: R₁, R₂ and inversion recovery
- `diffusion1d(filename)`: diffusion coefficient and hydrodynamic radius
- `tract(trosy, antitrosy)`: rotational correlation time from TRACT
- `calibration1d(filename)`: pulse-length calibration from a nutation experiment
- `kinetics1d(filename, times)`: intensities of named regions against time
- `r1rho([directory])`: R₁ρ relaxation dispersion
- `exchange1d([filenames])`: CEST and R₁ρ fitted by Bloch-McConnell simulation

# 2D experiments

Each opens a window to place and fit peaks, and returns the analysis when the window
closes: `results(expt)` gives one row per peak and `planeresults(expt)` one per peak per
plane. All take `skipplanes` to leave planes out of the fitting.

- `fit2d(files)`: positions, linewidths and amplitudes, with no model
- `relaxation2d(files)`: R₁ or R₂ from an exponential decay
- `recovery2d(files)`: inversion or saturation recovery
- `modelfit2d(files, x, equation, parameters)`: a model of your own
- `hetnoe2d(reference, saturated)`: heteronuclear NOE
- `ccr2d(decay, buildup, T)`: cross-correlated relaxation rates
- `methylccr2d(buildup, decay, T)`: methyl S²τc
- `cest2d(file)`: CEST profiles
- `cpmg2d(file)`: CPMG relaxation dispersion
- `peaktrack2d(files)`: peaks that move from plane to plane
- `titration2d(files)`: binding isotherms and a global Kd
- `rdc2d(; isotropic, aligned)`: one-bond couplings and RDCs
- `summaryplot(expt)`: a fitted parameter against residue number

# Other tools

- `viscosity(solvent, T)`: solvent viscosity
- `B1Calibration`: B₁ calibrations, from `calibration1d` or a reference pulse
"""
module NMRAnalysis

using LinearAlgebra: LAPACKException, PosDefException, SingularException
using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using REPL.TerminalMenus
using Reexport
using Statistics

include("fileselection.jl")
include("analyse.jl")
include("viscosity.jl")
include("output.jl")   # shared CSV column/value rules - see docs/src/advanced/conventions.md
include("fitting.jl")  # shared fitting utilities
include("b1.jl")       # B₁ calibration and inhomogeneity, shared by Analysis1D and Exchange1D
include("prompts.jl")  # parameter resolution, shared by Analysis1D and GUI2D

include("maybevector/MaybeVector.jl")
using .MaybeVectorModule

include("gui2d/GUI2D.jl")
using .GUI2D

include("analysis1d/Analysis1D.jl")
# N.B. deliberately no blanket `using .Analysis1D`: it exports `analyse`, which would
# collide with the registry-based `analyse` above. Names are brought in by the selective
# `@reexport using .Analysis1D: …` below instead.

using PrecompileTools
include("precompile.jl")

export analyse, register_analysis!, MultiFileRule
export viscosity
# B₁ fields: the calibration mapping power to field strength, and the distribution of that
# field across the sample (see src/b1.jl)
export B1Calibration, B1Distribution, b1average, b1rate, inhomogeneity, linearity
export shortpath

include("R1rho/R1rho.jl")
using .R1rho

include("exchange1d/Exchange1D.jl")
using .Exchange1D

@reexport using .MaybeVectorModule: MaybeVector, SingleElementVector, StandardVector
@reexport using .GUI2D: fit2d, relaxation2d, recovery2d, modelfit2d # IntensityExperiment
@reexport using .GUI2D: peaktrack2d, rdc2d, titration2d # MovingExperiment
@reexport using .GUI2D: hetnoe2d # HetNOEExperiment
@reexport using .GUI2D: cest2d # CESTExperiment
@reexport using .GUI2D: cpmg2d # CPMGExperiment
@reexport using .GUI2D: ccr2d # CCRExperiment
@reexport using .GUI2D: methylccr2d # methyl CCR (buildup/decay ratio)
@reexport using .GUI2D: summaryplot, results, planeresults

@reexport using .R1rho: r1rho, setupR1rhopowers

@reexport using .Exchange1D: exchange1d

# 1D analysis framework (Analysis1D) - the interactive replacements for the former
# readline-driven 1D routines. `analyse` is not re-exported (it would collide with the
# registry-based `analyse` above); use `analyse1d`.
@reexport using .Analysis1D: Region, Dataset1D, analyse1d, gui!, pickregion
@reexport using .Analysis1D: RelaxationExperiment, TractExperiment, NutationExperiment
@reexport using .Analysis1D: KineticsExperiment, DiffusionExperiment
@reexport using .Analysis1D: relaxation1d, tract, calibration1d, diffusion1d, kinetics1d
@reexport using .Analysis1D: calibrationanalysis
# results, and the series models a caller can choose between
@reexport using .Analysis1D: RegionResult, param
@reexport using .Analysis1D: SeriesModel, CurveFitModel, NoFitting, ExponentialModel

@info "NMRAnalysis.jl v$(pkgversion(NMRAnalysis)): type ?NMRAnalysis for a list of analyses"

end
