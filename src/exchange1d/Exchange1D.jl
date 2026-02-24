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
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using Plots
using REPL.TerminalMenus
using Statistics

# Import from parent module
using ..NMRAnalysis: analyse, register_analysis!, MultiFileRule
using ..NMRAnalysis: get1dregionandnoise, relaxation1d

# Include submodules in dependency order
include("fitting-with-errors.jl")
include("types.jl")
include("misc.jl")
include("models.jl")
include("experiments.jl")
include("liouvillian.jl")
include("params.jl")
include("problem.jl")

# """
#     exchange1d(filenames::Vector{String})

# Analyze 1D chemical exchange experiments (R1ρ and/or CEST).

# # Arguments
# - `filenames`: Vector of paths to NMR experiment directories

# # Returns
# Named tuple with:
# - `results`: Dict mapping models to FitResult
# - `plt`: Summary plot
# """
# function exchange1d(filenames::Vector{String})
#     @info "Starting exchange analysis..."

#     # 1. Classify experiments
#     r1rho_files = filter_r1rho_files(filenames)
#     cest_files = filter_cest_files(filenames)
#     r1_files = filter_r1_calibration_files(filenames)

#     @info "Found $(length(r1rho_files)) R1ρ, $(length(cest_files)) CEST, $(length(r1_files)) R1 calibration experiments"

#     if isempty(r1rho_files) && isempty(cest_files)
#         error("No R1ρ or CEST experiments found in provided files")
#     end

#     # 2. Load data
#     r1rho_data = isempty(r1rho_files) ? nothing : load_r1rho_data(r1rho_files)
#     return cest_data = isempty(cest_files) ? nothing : load_cest_data(cest_files)
# end

# # Registration with analysis system
# function __init__()
#     rule = MultiFileRule(expts -> begin
#                              oneD = filter(e -> "1d" in e.types, expts)
#                              r1rho = filter(e -> "r1rho" in e.types, oneD)
#                              cest = filter(e -> "cest" in e.types, oneD)
#                              r1cal = filter(e -> "relaxation" in e.types &&
#                                                 "R1" in e.features, oneD)
#                              combined = vcat(r1rho, cest, r1cal)
#                              # Only match if we have R1ρ or CEST experiments
#                              (length(r1rho) > 0 || length(cest) > 0) ? combined : nothing
#                          end,
#                          expts -> exchange1d([e.filename for e in expts]),
#                          "Exchange analysis (R1ρ and CEST)")
#     return register_analysis!(rule)
# end

end # module
