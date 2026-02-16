# Experiment types for exchange analysis

"""
    AbstractExperiment

Abstract type for NMR exchange experiments.
"""
abstract type AbstractExperiment end

"""
    R1rhoExperiment <: AbstractExperiment

R1ρ relaxation dispersion experiment data.

# Fields
- `Ω::Vector{Float64}`: Carrier offsets from resonance (rad/s)
- `ω1::Vector{Float64}`: Spin-lock field strengths (rad/s)
- `t::Vector{Float64}`: Relaxation times (s)
- `intensity::Vector{Float64}`: Observed intensities (empty for proposals)
- `σ::Vector{Float64}`: Measurement uncertainties (empty for proposals)
- `bf::Float64`: Spectrometer frequency (Hz)
"""
struct R1rhoExperiment <: AbstractExperiment
    Ω::Vector{Float64}
    ω1::Vector{Float64}
    t::Vector{Float64}
    intensity::Vector{Float64}
    σ::Vector{Float64}
    bf::Float64
end

"""
    R1rhoExperiment(Ω, ω1, t, bf)

Create an R1rhoExperiment without observed data (for simulation/proposal).
"""
function R1rhoExperiment(Ω::Vector{<:Real}, ω1::Vector{<:Real}, t::Vector{<:Real}, bf::Real)
    return R1rhoExperiment(Float64.(Ω), Float64.(ω1), Float64.(t), Float64[], Float64[],
                           Float64(bf))
end

"""
    CESTExperiment <: AbstractExperiment

Chemical Exchange Saturation Transfer experiment data.

# Fields
- `Ω::Vector{Float64}`: Saturation offsets from resonance (rad/s)
- `ω1::Vector{Float64}`: Saturation field strengths (rad/s)
- `t::Vector{Float64}`: Saturation times (s)
- `intensity::Vector{Float64}`: Observed intensities (empty for proposals)
- `σ::Vector{Float64}`: Measurement uncertainties (empty for proposals)
- `bf::Float64`: Spectrometer frequency (Hz)
"""
struct CESTExperiment <: AbstractExperiment
    Ω::Vector{Float64}
    ω1::Vector{Float64}
    t::Vector{Float64}
    intensity::Vector{Float64}
    σ::Vector{Float64}
    bf::Float64
end

"""
    CESTExperiment(Ω, ω1, t, bf)

Create a CESTExperiment without observed data (for simulation/proposal).
"""
function CESTExperiment(Ω::Vector{<:Real}, ω1::Vector{<:Real}, t::Vector{<:Real}, bf::Real)
    return CESTExperiment(Float64.(Ω), Float64.(ω1), Float64.(t), Float64[], Float64[],
                          Float64(bf))
end

"""
    hasdata(exp::AbstractExperiment)

Check if an experiment has observed data.
"""
hasdata(exp::AbstractExperiment) = !isempty(exp.intensity)

"""
    npoints(exp::AbstractExperiment)

Return the number of data points in the experiment.
"""
npoints(exp::AbstractExperiment) = length(exp.Ω)
