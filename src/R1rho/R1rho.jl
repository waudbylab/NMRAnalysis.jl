module R1rho

using CairoMakie
using Distributions: cdf, FDist
using GLMakie
using LinearAlgebra
using LsqFit
using Measurements: Measurements
using MonteCarloMeasurements: ±, pmean, pstd, register_primitive
using NMRTools
using Printf
using Random
using Statistics

export r1rho, setupR1rhopowers
using ..NMRAnalysis: select_expts, analyse_1d_calibration
using ..NMRAnalysis: register_analysis!, MultiFileRule

include("dataset.jl")
include("power.jl")
include("experiments.jl")
include("fitting.jl")
include("Kandkex.jl")
include("state.jl")
include("gui.jl")
include("calibrations.jl")

"""
    r1rho(directory_path=""; minvSL=250, maxvSL=1e6, scalefactor=:automatic)
    r1rho(filenames::Vector{String}; minvSL=250, maxvSL=1e6, scalefactor=:automatic)

Launch the R1rho analysis GUI for a given directory or experiment folder.

- `directory_path`: Path to the experiment folder or parent directory. If not specified, a dialog will prompt for selection.
- `filenames`: Vector of experiment files to analyze.
- `minvSL`: Minimum value for on-resonance spinlock strengths (in Hz). Default is 250 Hz. Set to zero to disable the lower limit.
  This allows low spinlock strengths with poor alignment to be excluded from the analysis.
- `maxvSL`: Maximum value for on-resonance spinlock strengths (in Hz). Default is 1 MHz (effectively no upper limit).
  This allows high spinlock strengths that may have detuning effects to be excluded from the analysis.
- `scalefactor`: Adjusts the size of the display window. By default, this is `2` for high-resolution displays and `1` for low-resolution displays. Omit or set to `:automatic` to use the default, or provide a numeric value to override.

# Example

```julia
# analyse experiments 11 and 12
r1rho(["examples/R1rho/11", "examples/R1rho/12"])

# open a dialog to select a directory, and scale the display size by 1.5
r1rho(scalefactor=1.5)

# filter out spinlock strengths below 500 Hz (e.g. to exclude poorly aligned experiments, 250 Hz is the default)
r1rho("examples/R1rho", minvSL=500)

# filter out spinlock strengths about 10 kHz (e.g. to exclude data with detuning effects)
r1rho("examples/R1rho", maxvSL=10_000)
```
"""
function r1rho(directory_path=""; minvSL=250, maxvSL=1e6, scalefactor=:automatic)
    filenames = select_expts(directory_path; experiment_type_filter="r1rho")
    isempty(filenames) && return
    return r1rho(filenames; minvSL=minvSL, maxvSL=maxvSL, scalefactor=scalefactor)
end

function r1rho(filenames::Vector{String}; minvSL=250, maxvSL=1e6, scalefactor=:automatic)
    if scalefactor == :automatic
        GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting")
    elseif scalefactor isa Number
        GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting",
                          scalefactor=scalefactor)
    else
        @error "scalefactor must be :auto or a number"
    end
    if minvSL >= maxvSL
        @error "minvSL ($minvSL Hz) must be less than maxvSL ($maxvSL Hz). Aborting R1rho analysis."
        return nothing
    end

    dataset = processexperiments(filenames; minvSL=minvSL, maxvSL=maxvSL)
    isnothing(dataset) && return nothing

    state = initialisestate(dataset)
    state[:filenames] = filenames
    return gui!(state)
end

function __init__()
    # Register primitives here - this gets called when the module is loaded
    # register the function with MonteCarloMeasurements as a primative
    register_primitive(safesqrt)

    # Register analysis rule for on-resonance R1rho experiments
    rule = MultiFileRule(expts -> begin
                             matched = filter(e -> "r1rho" in e.types &&
                                                       "1d" in e.types &&
                                                       "on_resonance" in
                                                       e.features,
                                              expts)
                             length(matched) > 0 ? matched : nothing
                         end,
                         expts -> r1rho([e.filename for e in expts]),
                         "On-resonance R1ρ dispersion")
    return register_analysis!(rule)
end

end