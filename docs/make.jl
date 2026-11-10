using Pkg
Pkg.activate(@__DIR__)  # Activate the docs environment
Pkg.instantiate()       # Install all dependencies

using Documenter, NMRAnalysis
using Markdown

ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

# Screenshots are made as the pages build, by opening the example analyses with
# NMRAnalysis.screenshot, so they always show the current windows. The pages' example blocks
# run in modules of their own, and reach these helpers as Main.<name>.

"Save a figure as a PNG beside the page."
savepng(file, fig) = NMRAnalysis.GLMakie.save(file, fig; backend=NMRAnalysis.GLMakie,
                                              px_per_unit=1)

"The first `n` rows of `results(expt)`, as a table."
function resultstable(rows; n=6)
    cell(x) = x isa AbstractFloat ? string(round(x; sigdigits=4)) : string(x)
    cols = keys(first(rows))
    lines = ["| " * join(cols, " | ") * " |", "|" * repeat("---|", length(cols))]
    for row in first(rows, n)
        push!(lines, "| " * join(map(cell, values(row)), " | ") * " |")
    end
    length(rows) > n && push!(lines, "| … |" * repeat(" |", length(cols) - 1))
    return Markdown.parse(join(lines, "\n"))
end

# DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRAnalysis); recursive=true)

makedocs(;
         modules=[NMRAnalysis],
         # the API pages hold every docstring in their module, so they outgrow the limit
         format=Documenter.HTML(; canonical="https://waudbylab.org/NMRAnalysis.jl/stable/",
                                size_threshold_ignore=["api/gui2d.md", "api/analysis1d.md"]),
         pages=["Home" => "index.md",
                "Quick start" => "quickstart.md",
                "Automatic analysis" => "analyses/analyse.md",
                "1D experiments" => ["The analysis window" => "analyses/1d/overview.md",
                                     "Relaxation (R₁, R₂)" => "analyses/1d/relaxation.md",
                                     "Diffusion" => "analyses/1d/diffusion.md",
                                     "TRACT" => "analyses/1d/tract.md",
                                     "Pulse Calibration" => "analyses/1d/calibration.md",
                                     "Kinetics" => "analyses/1d/kinetics.md",
                                     "Peak Tracking" => "analyses/1d/peaktracking.md"],
                "1D chemical exchange" =>
                    ["R1ρ Relaxation Dispersion" => "analyses/r1rho.md",
                     "Exchange Fitting (CEST, R1ρ, R1)" =>
                         ["Overview" => "analyses/exchange1d/overview.md",
                          "Exchange Models" => "analyses/exchange1d/models.md",
                          "Theory and Calculation Methods" => "analyses/exchange1d/theory.md"]],
                "2D experiments" => ["The analysis window" => "analyses/2d/overview.md",
                                     "Relaxation (R₁, R₂)" => "analyses/2d/relaxation.md",
                                     "Heteronuclear NOE" => "analyses/2d/hetnoe.md",
                                     "Titrations" => "analyses/2d/titration.md",
                                     "RDCs and J-Coupling" => "analyses/2d/rdc.md",
                                     "Magnetisation Recovery" => "analyses/2d/magnetisationrecovery.md",
                                     "Cross-Correlated Relaxation (CCR)" => "analyses/2d/ccr.md",
                                     "Methyl CCR (S²τc)" => "analyses/2d/methylccr.md",
                                     "Under Development" => ["CEST" => "analyses/2d/cest.md",
                                                             "CPMG Dispersion" => "analyses/2d/cpmg.md",
                                                             "Paramagnetic Relaxation Enhancement" => "analyses/2d/pre.md"],
                                     "Reference" =>
                                         ["Peak List Formats" => "analyses/2d/peaklistformats.md",
                                          "Summary Plots" => "analyses/2d/summary.md",
                                          "Simple Fitting" => "analyses/2d/fit.md",
                                          "Peak Tracking" => "analyses/2d/peaktracking.md",
                                          "Custom Models" => "analyses/2d/modelfit.md"]],
                "Tutorials" => ["¹⁹F R1ρ Acquisition & Analysis" => "tutorials/r1rho.md"],
                "Ecosystem" => "ecosystem.md",
                "Advanced" => ["How an Analysis Works" => "advanced/pipeline.md",
                               "Output and Interface Conventions" => "advanced/conventions.md",
                               "Analysis Rules" => "advanced/analysis_rules.md",
                               "Creating New 1D Analyses" => "advanced/creating_1d_analyses.md",
                               "Creating New 2D Analyses" => "advanced/creating_2d_analyses.md",
                               "Extending Exchange1D" => "advanced/extending_exchange1d.md",
                               "API" => ["NMRAnalysis" => "api.md",
                                         "Analysis1D" => "api/analysis1d.md",
                                         "GUI2D" => "api/gui2d.md",
                                         "R1ρ" => "api/r1rho.md",
                                         "Exchange1D" => "api/exchange1d.md",
                                         "PeakTrack1D" => "api/peaktrack1d.md",
                                         "MaybeVector" => "api/maybevector.md"],
                               "Index" => "indexes.md"]],
         sitename="NMRAnalysis.jl",
         authors="Chris Waudby",
         warnonly=[:missing_docs])

deploydocs(;
           repo="github.com/waudbylab/NMRAnalysis.jl.git",
           devbranch="main",
           # Build a browsable preview for every pull request, at
           # https://waudbylab.github.io/NMRAnalysis.jl/previews/PR<number>/
           push_preview=true)
