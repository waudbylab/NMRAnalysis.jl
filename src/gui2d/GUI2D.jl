module GUI2D

using CairoMakie
using DelimitedFiles
using GLMakie
using Graphs
using LinearAlgebra: ColumnNorm, Symmetric, cholesky, diag, issuccess, qr
using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using OrderedCollections
using PrettyTables
using REPL.TerminalMenus
using Statistics
using ..MaybeVectorModule
# shared output rules - see src/output.jl and docs/src/advanced/conventions.md
using ..NMRAnalysis: stderrors   # standard errors that survive a singular covariance
using ..NMRAnalysis: csvcolumn, csvcolumns, csvvalue, safename, sanitizelabel, backupfile,
                     backupfolder, writetable
# parameter resolution: argument, then annotation/acqus, then ask - see src/prompts.jl
using ..NMRAnalysis: annotation, acqusvalue, ask, askvector
using ..NMRAnalysis: B1Calibration
# the call that produced an analysis, for summary.txt - see src/calls.jl
using ..NMRAnalysis: AnalysisCall, analysiscall, callstring

include("util.jl")
include("types.jl")
include("parameters.jl")
include("specdata.jl")
include("peaks.jl")
include("experiments.jl")
include("fitting.jl")
include("models.jl")
include("clustering.jl")
include("state.jl")
include("gui.jl")
include("mouse.jl")
include("keyboard.jl")
include("files.jl")
include("output.jl")
include("visualisation.jl")
include("summary.jl")

export MaybeVector, SingleElementVector, StandardVector
# export gui!

# export IntensityExperiment
export fit2d
export relaxation2d
export recovery2d
export modelfit2d

# export MovingExperiment
export peaktrack2d
export rdc2d
export titration2d

# export HetNOEExperiment
export hetnoe2d

# export CESTExperiment
export cest2d

# export CPMGExperiment
export cpmg2d

# export CCRExperiment
export ccr2d

# methyl CCR (buildup/decay ratio, eq 7)
export methylccr2d

# results-summary plotting
export summaryplot

# what an analysis returns when its window closes
export results, planeresults

end
