"""
    Analysis1D

Unified framework for lightweight 1D NMR analyses (relaxation, TRACT, nutation
calibration, STD, kinetics, …). See `PLAN.md` in this directory for the design.

The analysis core operates on plain `Trace`/`Planes`/`Dataset1D` values and has no GUI
or NMRData dependency in its computational path (the agreed "keep the science pure"
split). NMRData is touched only by the loader adapters in `loaders.jl`. The interactive
GUI is a later phase layered on top.
"""
module Analysis1D

using CairoMakie
using GLMakie
using LsqFit
using Measurements
using NMRTools
using Statistics

# pure analysis core (no Makie dependency in these files)
include("types.jl")
include("reductions.jl")
include("seriesmodels.jl")
include("experiments.jl")
include("std.jl")
include("loaders.jl")

# interactive GUI
include("visualisation.jl")
include("state.jl")
include("gui.jl")

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

# interactive GUI
export gui!

# file loaders (stdnmr rather than std, to avoid colliding with Statistics.std)
export relaxation, tract, nutation, stdnmr, kinetics

end # module Analysis1D
