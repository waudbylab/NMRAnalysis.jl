module GUI2D

using CairoMakie
using DelimitedFiles
using GLMakie
using Graphs
using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using OrderedCollections
using PrettyTables
using REPL.TerminalMenus
using Statistics
using ..MaybeVectorModule

include("util.jl")
include("types.jl")
include("parameters.jl")
include("specdata.jl")
include("peaks.jl")
include("experiments.jl")
include("models.jl")
include("clustering.jl")
include("state.jl")
include("gui.jl")
include("mouse.jl")
include("keyboard.jl")
include("files.jl")
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

end