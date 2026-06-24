"""
    Analysis1D

Unified framework for lightweight 1D NMR analyses (relaxation, TRACT, nutation
calibration, STD, kinetics, …). See `PLAN.md` in this directory for the design.

The analysis core operates on plain `Trace`/`Planes`/`Dataset1D` values and has no GUI
or NMRData dependency in its computational path, so it is testable headless. NMRData is
touched only by the `load_*` adapters in `loaders.jl`. The interactive GUI is a later
phase layered on top.
"""
module Analysis1D

using LsqFit
using Measurements
using NMRTools
using Statistics

include("types.jl")
include("reductions.jl")
include("seriesmodels.jl")
include("experiments.jl")
include("std.jl")
include("loaders.jl")

"""
    analyse1d(experiment) -> NamedTuple

Run a 1D analysis. Equivalent to [`analyse`](@ref) on an `Experiment1D`; provided under
a distinct name so it can be re-exported without colliding with the registry-based
`analyse` dispatcher in the parent module.
"""
analyse1d(e) = analyse(e)

# core data types
export Trace, Planes, Region, Dataset1D
export column, hasvar, nplanes, groupseries

# reductions & models
export Reduction, Integrate, integrate, integrals
export SeriesModel, CurveFitModel, NoFitting, ContrastModel
export ExponentialModel, RecoveryModel, DampedSinusoidModel

# experiments
export Experiment1D, analyse, analyse1d
export RelaxationExperiment, TractExperiment, NutationExperiment, KineticsExperiment
export STDExperiment
export SeriesResult, param

# file loaders
export load_relaxation, load_tract, load_nutation, load_std, load_kinetics

end # module Analysis1D
