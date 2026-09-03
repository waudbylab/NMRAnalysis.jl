using Pkg
Pkg.activate(@__DIR__)  # Activate the docs environment
Pkg.instantiate()       # Install all dependencies

using Documenter, NMRAnalysis

ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

# DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRAnalysis); recursive=true)

makedocs(;
         modules=[NMRAnalysis],
         format=Documenter.HTML(; canonical="https://waudbylab.org/NMRAnalysis.jl/stable/"),
         pages=["Home" => "index.md",
                "Quick Start" => "quickstart.md",
                "Automatic Analysis" => "analyses/analyse.md",
                "1D Experiments" => ["Diffusion" => "analyses/diffusion.md",
                                     "Relaxation" => "analyses/relaxation.md",
                                     "TRACT" => "analyses/tract.md"],
                "1D Chemical Exchange" => ["Calibration" => "analyses/calibration.md",
                                           "R1ρ Relaxation Dispersion" => "analyses/r1rho.md",
                                           "Exchange Fitting (CEST, R1ρ, R1)" => "analyses/exchange1d.md"],
                "2D Experiments" => ["Overview" => "analyses/2d/overview.md",
                                     "Relaxation (R₁, R₂)" => "analyses/2d/relaxation.md",
                                     "Heteronuclear NOE" => "analyses/2d/hetnoe.md",
                                     "Titrations" => "analyses/2d/titration.md",
                                     "RDCs and J-Coupling" => "analyses/2d/rdc.md",
                                     "Magnetisation Recovery" => "analyses/2d/magnetisationrecovery.md",
                                     "Cross-Correlated Relaxation (CCR)" => "analyses/2d/ccr.md",
                                     "Methyl CCR (S²τc)" => "analyses/2d/methylccr.md",
                                     "Under Development" => ["CEST" => "analyses/2d/cest.md",
                                                             "CPMG Dispersion" => "analyses/2d/cpmg.md"],
                                     "Reference" => ["Peak List Formats" => "analyses/2d/peaklistformats.md",
                                                     "Summary Plots" => "analyses/2d/summary.md",
                                                     "Simple Fitting" => "analyses/2d/fit.md",
                                                     "Peak Tracking" => "analyses/2d/peaktracking.md",
                                                     "Custom Models" => "analyses/2d/modelfit.md"]],
                "Tutorials" => ["¹⁹F R1ρ Acquisition & Analysis" => "tutorials/r1rho.md"],
                "Ecosystem" => "ecosystem.md",
                "Advanced" => ["Analysis Rules" => "advanced/analysis_rules.md",
                               "Creating New 1D Analyses" => "advanced/creating_1d_analyses.md",
                               "Creating New 2D Analyses" => "advanced/creating_2d_analyses.md",
                               "Extending Exchange1D" => "advanced/extending_exchange1d.md",
                               "API" => ["NMRAnalysis" => "api.md",
                                         "GUI2D" => "api/gui2d.md",
                                         "R1ρ" => "api/r1rho.md",
                                         "Exchange1D" => "api/exchange1d.md",
                                         "MaybeVector" => "api/maybevector.md"],
                               "Index" => "indexes.md"]],
         sitename="NMRAnalysis.jl",
         authors="Chris Waudby",
         warnonly=[:missing_docs],)

deploydocs(;
           repo="github.com/waudbylab/NMRAnalysis.jl.git",
           devbranch="main",)
