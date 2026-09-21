"""
    Analysis1D

Unified framework for lightweight 1D NMR analyses (relaxation, TRACT, nutation
calibration, kinetics, …). See `PLAN.md` in this directory for the design.

The analysis core operates on plain `Trace`/`Planes`/`Dataset1D` values and has no GUI
or NMRData dependency in its computational path (the agreed "keep the science pure"
split). NMRData is touched only by the adapters in `nmrdata.jl`.

Each analysis lives in a single self-contained `expt-<name>.jl` holding its entry point,
type, interface methods, science and presentation - the layout used by `GUI2D` and
`Exchange1D`. `experiments.jl` documents the interface and includes them; see
`expt-tract.jl` for the fullest example.
"""
module Analysis1D

using CairoMakie
using GLMakie
using LsqFit
using Measurements
using NMRTools
using OrderedCollections
using REPL.TerminalMenus
using Statistics

using ..NMRAnalysis: register_analysis!, MultiFileRule, viscosity
using ..NMRAnalysis: stderrors   # standard errors that survive a singular covariance
# B₁ calibration: the type and its accessors, built from a nutation fit in expt-nutation.jl.
# `import` rather than `using`, because expt-nutation.jl adds constructor methods to it:
# `using` brings the name in to call, not to extend.
import ..NMRAnalysis: B1Calibration
using ..NMRAnalysis: refpower, ν1ref, linearity, inhomogeneity, powers, fields
# shared output rules - see src/output.jl and docs/src/advanced/conventions.md
using ..NMRAnalysis: csvcolumn, csvcolumns, csvvalue, safename, sanitizelabel, backupfile,
                     backupfolder, writetable, shortpath

# the analysis core: no Makie in the computational path. (An experiment that saves a
# figure of its own draws it in its `expt-*.jl`, as Exchange1D's experiments do, so that
# everything about one analysis stays in one file.)
include("types.jl")
include("integration.jl")
include("seriesmodels.jl")
include("nmrdata.jl")
include("prompts.jl")   # parameter resolution: argument, then annotation/acqus, then ask
include("experiments.jl")   # interface + pipeline; includes one expt-*.jl per experiment
include("files.jl")

# interactive GUI
include("visualisation.jl")
include("state.jl")
include("gui.jl")

"""
    analyse1d(experiment) -> Vector{RegionResult}

Run a 1D analysis. Equivalent to [`analyse`](@ref) on an `Experiment1D`; provided under
a distinct name so it can be re-exported without colliding with the registry-based
`analyse` dispatcher in the parent module.
"""
analyse1d(e) = analyse(e)

# core data types
export Trace, Planes, Region, Dataset1D
export column, hasvar, nplanes, groupseries

# measurement & models
export integrate
export SeriesModel, CurveFitModel, NoFitting
export ExponentialModel, RecoveryModel, DampedSinusoidModel, StejskalTannerModel

# experiments
export Experiment1D, analyse, analyse1d, run1d, Integration
export RelaxationExperiment, TractExperiment, NutationExperiment, KineticsExperiment
export DiffusionExperiment
export RegionResult, SeriesResult, param

# interactive GUI
export gui!, pickregion

# top-level entry points (each opens the GUI, or analyses directly given an
# `integration` triple).
export relaxation1d, tract, calibration1d, diffusion1d, kinetics1d
# fit a calibration without a window, and save it where the caller decides
export calibrationanalysis

"""
Register the interactive 1D analyses with the analysis-dispatch registry, so `analyse`
routes annotated experiments here. These replace the registrations previously made by the
readline-based routines.

Only the analyses that can run from a filename alone are registered: TRACT needs a
TROSY/anti-TROSY pair, diffusion needs the gradient list, and kinetics needs named
regions, so those are invoked directly rather than dispatched.
"""
function __init__()
    register_analysis!(["1d", "relaxation"], ["R1"],
                       e -> relaxation1d(e.filename), "1D R1 relaxation")
    register_analysis!(["1d", "relaxation"], ["R2"],
                       e -> relaxation1d(e.filename), "1D R2 relaxation")
    # Every nutation calibration selected at once, since several power levels make a
    # calibration curve where one makes only a point on it.
    isnutation(e) = "1d" in e.types && "calibration" in e.types && "nutation" in e.features
    function nutations(expts)
        matched = filter(isnutation, expts)
        return isempty(matched) ? nothing : matched
    end
    return register_analysis!(MultiFileRule(nutations,
                                            expts -> calibration1d([e.filename
                                                                    for e in expts]),
                                            "1D nutation calibration"))
end

end # module Analysis1D
