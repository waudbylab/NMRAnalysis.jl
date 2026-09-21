"""
    Exchange1D

Module for 1D chemical exchange analysis using Bloch-McConnell equations.

Supports R1ρ relaxation dispersion and CEST experiments with:
- NoExchange (null model), TwoState, and TwoStateBinding exchange models
- Full Bloch-McConnell Liouvillian construction
- Joint fitting of multiple experiment types
"""
module Exchange1D

using CairoMakie
using ComponentArrays
using GLMakie
using InteractiveUtils: subtypes
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using PrettyTables
using REPL.TerminalMenus
using Statistics

# Import from parent module
import ..NMRAnalysis  # module itself, for pkgversion(NMRAnalysis)
using ..NMRAnalysis: analyse, register_analysis!, MultiFileRule
using ..NMRAnalysis: select_expts
using ..NMRAnalysis: stderrors   # standard errors that survive a singular covariance
# B₁ field strength and its distribution across the sample - see src/b1.jl
using ..NMRAnalysis: B1Calibration, B1Distribution, b1average, b1rate, inhomogeneity,
                     npoints
# shared output rules - see src/output.jl and docs/src/advanced/conventions.md
using ..NMRAnalysis: csvcolumn, csvvalue, safename, backupfile, backupfolder,
                     writetable
# Interactive region selection for `integrate!` (replaces the former readline prompts).
using ..Analysis1D: pickregion, calibrationanalysis

# Include submodules in dependency order
include("fitting-with-errors.jl")
include("types.jl")
include("misc.jl")
include("models.jl")
include("experiments.jl")
include("liouvillian.jl")
include("params.jl")
include("problem.jl")
include("plots.jl")
include("overlay.jl")
include("interface.jl")
include("results.jl")
include("files.jl")

export exchange1d

# Registration with analysis system
"""A nutation calibration among the selected experiments is not something to fit: it is
where the spin-lock field strengths and the B₁ inhomogeneity come from. Matched alongside
the experiments so that a folder holding both is analysed with measured fields rather than
nominal ones."""
iscalibration(e) = "calibration" in e.types && "nutation" in e.features

"""The experiments an exchange analysis can be run on, or `nothing` where the selection
holds nothing to fit. A CEST or off-resonance R1ρ experiment is what makes an exchange
analysis worth offering; R1 and on-resonance R1ρ experiments join a fit but do not trigger
one, and a calibration only ever joins."""
function exchangeexperiments(expts)
    oneD = filter(e -> "1d" in e.types, expts)
    cest = filter(e -> "cest" in e.types, oneD)
    r1cal = filter(e -> "relaxation" in e.types && "R1" in e.features, oneD)
    onres = filter(e -> "r1rho" in e.types && "on_resonance" in e.features, oneD)
    offres = filter(e -> "r1rho" in e.types && "off_resonance" in e.features, oneD)
    isempty(cest) && isempty(offres) && return nothing
    return vcat(cest, r1cal, onres, offres, filter(iscalibration, oneD))
end

"""Run the analysis the dispatcher matched, with any calibration experiments among the
selection used as the calibration rather than fitted as data."""
function runexchange(expts)
    calibration = [e.filename for e in expts if iscalibration(e)]
    return exchange1d([e.filename for e in expts if !iscalibration(e)];
                      calibration=isempty(calibration) ? nothing : calibration)
end

function __init__()
    return register_analysis!(MultiFileRule(exchangeexperiments, runexchange,
                                            "Exchange analysis (CEST / R1rho)"))
end

end # module
