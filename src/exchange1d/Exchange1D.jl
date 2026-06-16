"""
    Exchange1D

Module for 1D chemical exchange analysis using Bloch-McConnell equations.

Supports R1ρ relaxation dispersion and CEST experiments with:
- NoExchange (null model), TwoState, and TwoStateBinding exchange models
- Full Bloch-McConnell Liouvillian construction
- Joint fitting of multiple experiment types
"""
module Exchange1D

using ComponentArrays
using InteractiveUtils: subtypes
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using Plots
using PrettyTables
using REPL.TerminalMenus
using Statistics

# Import from parent module
using ..NMRAnalysis: analyse, register_analysis!, MultiFileRule
using ..NMRAnalysis: get1dregionandnoise, relaxation1d
using ..NMRAnalysis: select_expts

# Include submodules in dependency order
include("fitting-with-errors.jl")
include("types.jl")
include("misc.jl")
include("models.jl")
include("experiments.jl")
include("liouvillian.jl")
include("params.jl")
include("problem.jl")
include("interface.jl")
include("results.jl")

# Registration with analysis system
function __init__()
    rule = MultiFileRule(expts -> begin
                             oneD = filter(e -> "1d" in e.types, expts)
                             cest = filter(e -> "cest" in e.types, oneD)
                             r1cal = filter(e -> "relaxation" in e.types &&
                                                "R1" in e.features, oneD)
                             onres = filter(e -> "r1rho" in e.types &&
                                                "on_resonance" in e.features, oneD)
                             offres = filter(e -> "r1rho" in e.types &&
                                                "off_resonance" in e.features, oneD)
                             combined = vcat(cest, r1cal, onres, offres)
                             (length(cest) > 0 || length(offres) > 0) ? combined : nothing
                         end,
                         expts -> exchange1d([e.filename for e in expts]),
                         "Exchange analysis (CEST / R1rho)")
    return register_analysis!(rule)
end

end # module
