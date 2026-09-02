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
include("diffusion1d.jl")
include("relaxation1d.jl")
include("tract.jl")

include("maybevector/MaybeVector.jl")
using .MaybeVectorModule

include("gui2d/GUI2D.jl")
using .GUI2D

include("analysis1d/Analysis1D.jl")
# N.B. deliberately no blanket `using .Analysis1D`: it exports `analyse`, which would
# collide with the registry-based `analyse` above. Names are brought in by the selective
# `@reexport using .Analysis1D: …` below instead.

using PrecompileTools
include("precompile.jl")

export analyse, register_analysis!, MultiFileRule
export viscosity
export diffusion1d
export relaxation1d
export tract

include("R1rho/R1rho.jl")
using .R1rho

include("exchange1d/Exchange1D.jl")
using .Exchange1D

@reexport using .MaybeVectorModule: MaybeVector, SingleElementVector, StandardVector
@reexport using .GUI2D: fit2d, relaxation2d, recovery2d, modelfit2d # IntensityExperiment
@reexport using .GUI2D: peaktrack2d, rdc2d, titration2d # MovingExperiment
@reexport using .GUI2D: hetnoe2d # HetNOEExperiment
@reexport using .GUI2D: cest2d # CESTExperiment
@reexport using .GUI2D: cpmg2d # CPMGExperiment
@reexport using .GUI2D: ccr2d # CCRExperiment
@reexport using .GUI2D: methylccr2d # methyl CCR (buildup/decay ratio)
@reexport using .GUI2D: summaryplot

@reexport using .R1rho: r1rho, setupR1rhopowers

@reexport using .Exchange1D: exchange1d

# 1D analysis framework (Analysis1D). `analyse` is not re-exported (would collide with the
# registry-based `analyse` above; use `analyse1d`). `tract` is not re-exported either, as it
# collides with the legacy top-level `tract`; reach the new one as `Analysis1D.tract`.
@reexport using .Analysis1D: Region, Dataset1D, analyse1d, gui!
@reexport using .Analysis1D: RelaxationExperiment, TractExperiment, NutationExperiment
@reexport using .Analysis1D: KineticsExperiment, STDExperiment
@reexport using .Analysis1D: relaxation, nutation, stdnmr, kinetics

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
- diffusion1d([filename])
- tract([trosy_filename, antitrosy_filename])
- r1rho([directory_path]; minvSL=250, maxvSL=1e6, scalefactor=:automatic)
- exchange1d([filenames]) - CEST / R1ρ chemical exchange analysis

# 1D Analysis Framework (interactive; build an experiment then `gui!(expt)`)

- relaxation(spec; ir=false)        |> gui!
- nutation(spec)                    |> gui!
- stdnmr(spec, sat, tsat; regions)  |> gui!
- kinetics(spec, times; regions)    |> gui!
- Analysis1D.tract(trosy, antitrosy) |> gui!

# 2D Experiment Analysis Routines

- fit2d(inputfilenames)
- relaxation2d(inputfilenames, relaxationtimes | taufilename)
- recovery2d(inputfilenames, relaxationtimes | taufilename)
- modelfit2d(inputfilenames, xvalues, equation, parameters)
- hetnoe2d(inputfilenames, saturationlist)
- ccr2d(decay_experiments, buildup_experiments, Trelax)
- methylccr2d(buildup_experiment, decay_experiment, T; C=3/4)
- cest2d(inputfilenames; B1, Tsat)
- cpmg2d(inputfilename; Trelax, vCPMG | ncyc)

Current working directory: $(pwd())
"""

end
