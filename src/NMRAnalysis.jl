module NMRAnalysis

using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using Plots
using REPL.TerminalMenus
using Reexport
using Statistics

include("fileselection.jl")
include("regions1d.jl")
include("analyse.jl")
include("viscosity.jl")
include("calibration1d.jl")
include("diffusion.jl")
include("relaxation1d.jl")
include("tract.jl")

include("gui2d/GUI2D.jl")
using .GUI2D

using PrecompileTools
include("precompile.jl")

export analyse, register_analysis!, MultiFileRule
export viscosity
export diffusion
export relaxation1d
export tract

include("R1rho/R1rho.jl")
using .R1rho

@reexport using .GUI2D: MaybeVector
@reexport using .GUI2D: intensities2d, relaxation2d, recovery2d, modelfit2d # IntensityExperiment
@reexport using .GUI2D: hetnoe2d # HetNOEExperiment
@reexport using .GUI2D: cest2d # CESTExperiment
@reexport using .GUI2D: cpmg2d # CPMGExperiment
@reexport using .GUI2D: pre2d # PREExperiment

@reexport using .R1rho: r1rho, setupR1rhopowers

@info """
NMRAnalysis.jl (v$(pkgversion(NMRAnalysis)))

1. set your working directory to a convenient location, e.g.
   cd("/Users/chris/NMR/crick-702/my_experiment_directory")
2. call the desired analysis routine
3. use `?function_name` to get help on any function

# Generic Analysis (alpha)

- analyse(filename)
- analyse([filename1, filename2, ...])

# 1D Experiment Analysis Routines

- relaxation1d([filename])
- diffusion([filename])
- tract([trosy_filename, antitrosy_filename])
- r1rho([directory_path]; minvSL=250, maxvSL=1e6, scalefactor=:automatic)

# 2D Experiment Analysis Routines

- intensities2d(inputfilenames)
- relaxation2d(inputfilenames, relaxationtimes | taufilename)
- recovery2d(inputfilenames, relaxationtimes | taufilename)
- modelfit2d(inputfilenames, xvalues, equation, parameters)
- hetnoe2d(inputfilenames, saturationlist)
- cest2d(inputfilenames; B1, Tsat)
- cpmg2d(inputfilename; Trelax, vCPMG | ncyc)

Current working directory: $(pwd())
"""

end
