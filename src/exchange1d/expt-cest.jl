"""
    CESTExperiment <: AbstractExperiment

Chemical Exchange Saturation Transfer experiment. Stores saturation offsets,
B1 field strength, saturation time, and observed/predicted intensity profiles.

Loading, integration, and simulation are not yet implemented.
"""
struct CESTExperiment <: AbstractExperiment
    spec::Any                # NMRData (untyped for testability)
    field_teslas::Float64
    sampleconcentrations::Dict{String,Float64}
    offsets::Vector{Float64}                    # saturation offsets in Hz
    B1_field::Float64                           # saturation field strength in Hz
    saturation_time::Float64                    # seconds
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
end

"""
    default_spin_params(expt::CESTExperiment, N::Int) -> Vector{Pair{Symbol,Any}}

Return spin parameter entries needed by this CEST experiment: R1, R2 for the
experiment's field, and chemical shifts (delta).
"""
function default_spin_params(expt::CESTExperiment, N::Int)
    fl = field_label(expt)
    return Pair{Symbol,Any}[:delta => zeros(N),
                            Symbol("R2_", fl) => fill(10.0, N),
                            Symbol("R1_", fl) => [1.5]]
end

"""
    default_nuisance_params(expt::CESTExperiment) -> Vector{Pair{Symbol,Any}}

Return flat nuisance parameter entries for this CEST experiment, tagged with
the experiment type and field, e.g. `:CEST_14p1T_I0`.
"""
function default_nuisance_params(expt::CESTExperiment)
    tag = Symbol("CEST_", field_label(expt))
    return Pair{Symbol,Any}[Symbol(tag, "_I0") => 1.0]
end
