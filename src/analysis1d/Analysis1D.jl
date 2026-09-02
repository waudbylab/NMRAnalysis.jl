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

using ..NMRAnalysis: register_analysis!, viscosity

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

export StejskalTannerModel

# experiments
export Experiment1D, analyse, analyse1d, run1d, Integration
export RelaxationExperiment, TractExperiment, NutationExperiment, KineticsExperiment
export DiffusionExperiment, STDExperiment
export SeriesResult, param

# interactive GUI
export gui!, pickregion

# top-level entry points (each opens the GUI, or analyses directly given an
# `integration` triple). `std1d` rather than `std` also avoids colliding with Statistics.
export relaxation1d, tract, calibration1d, diffusion1d, std1d, kinetics1d

"""
Register the interactive 1D analyses with the analysis-dispatch registry, so `analyse`
routes annotated experiments here. These replace the registrations previously made by the
readline-based routines.

Only the analyses that can run from a filename alone are registered: TRACT needs a
TROSY/anti-TROSY pair, diffusion needs the gradient list, and STD/kinetics need named
regions, so those are invoked directly rather than dispatched.
"""
function __init__()
    register_analysis!(["1d", "relaxation"], ["R1"],
                       e -> relaxation1d(e.filename), "1D R1 relaxation")
    register_analysis!(["1d", "relaxation"], ["R2"],
                       e -> relaxation1d(e.filename), "1D R2 relaxation")
    return register_analysis!(["1d", "calibration"], ["nutation"],
                              e -> calibration1d(e.filename), "1D nutation calibration")
end

end # module Analysis1D
