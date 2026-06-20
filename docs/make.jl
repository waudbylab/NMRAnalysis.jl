using Pkg
Pkg.activate(@__DIR__)  # Activate the docs environment
Pkg.instantiate()       # Install all dependencies

using Documenter, NMRAnalysis

ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

# DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRAnalysis); recursive=true)

makedocs(;
         modules=[NMRAnalysis],
         format=Documenter.HTML(),
         pages=["Home" => "index.md",
                "Quick Start" => "quickstart.md",
                "Analyses" => ["Diffusion (1D)" => "analyses/diffusion.md",
                               "Relaxation (1D)" => "analyses/relaxation.md",
                               "TRACT (1D)" => "analyses/tract.md",
                               "Calibration (1D)" => "analyses/calibration.md",
                               "R1ρ (1D)" => "analyses/r1rho.md",
                               "2D Experiments" => ["Overview" => "analyses/2d/overview.md",
                                   "Peak Fitting" => "analyses/2d/fit.md",
                                   "Relaxation (T1/T2)" => "analyses/2d/relaxation.md",
                                   "Recovery" => "analyses/2d/recovery.md",
                                   "Custom Model" => "analyses/2d/modelfit.md",
                                   "Heteronuclear NOE" => "analyses/2d/hetnoe.md",
                                   "CPMG Dispersion" => "analyses/2d/cpmg.md",
                                   "CEST" => "analyses/2d/cest.md",
                                   "PRE" => "analyses/2d/pre.md",
                                   "CCR" => "analyses/2d/ccr.md"]],
                "Tutorials" => ["19F R1ρ acquisition & analysis" => "tutorials/r1rho.md"],
                "Advanced" => ["Analysis Rules" => "advanced/analysis_rules.md",
                               "Creating new 2D analyses" => "advanced/creating_2d_analyses.md",
                               "API" => "api.md",
                               "Index" => "indexes.md"]],
         sitename="NMRAnalysis.jl",
         authors="Chris Waudby",
         warnonly=[:missing_docs],)

deploydocs(;
           repo="github.com/waudbylab/NMRAnalysis.jl.git",
           devbranch="main",)
