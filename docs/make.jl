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
                "1D Analysis" => ["Diffusion" => "analyses/diffusion.md",
                                  "Relaxation" => "analyses/relaxation.md",
                                  "TRACT" => "analyses/tract.md",
                                  "Calibration" => "analyses/calibration.md",
                                  "R1ρ" => "analyses/r1rho.md"],
                "2D Analysis" => ["Overview" => "analyses/2d/overview.md",
                                  "Peak List Formats" => "analyses/2d/peaklistformats.md",
                                  "Summary Plots" => "analyses/2d/summary.md",
                                  "Simple Fitting" => "analyses/2d/fit.md",
                                  "Peak Tracking" => "analyses/2d/peaktracking.md",
                                  "CEST" => "analyses/2d/cest.md",
                                  "CPMG Dispersion" => "analyses/2d/cpmg.md",
                                  "Cross-Correlated Relaxation (CCR)" => "analyses/2d/ccr.md",
                                  "Heteronuclear NOE" => "analyses/2d/hetnoe.md",
                                  "Magnetisation Recovery" => "analyses/2d/magnetisationrecovery.md",
                                  "Methyl CCR (S²τc)" => "analyses/2d/methylccr.md",
                                  "PREs" => "analyses/2d/pre.md",
                                  "RDCs and J-Coupling" => "analyses/2d/rdc.md"],
                "Relaxation (R₁, R₂)" => "analyses/2d/relaxation.md",
                "Titrations" => "analyses/2d/titration.md",
                "Custom Models" => "analyses/2d/modelfit.md",
                "Tutorials" => ["¹⁹F R1ρ acquisition & analysis" => "tutorials/r1rho.md"],
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
