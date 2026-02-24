function componentnames(spec::NMRData)
    components = sample(spec, :sample, :components)
    return [component["name"] for component in components]
end
componentnames(expt::AbstractExperiment) = componentnames(expt.spec)

"""
    sampleconcentrations(spec::NMRData) -> Dict{String,Float64}

Extract a Dict mapping molecule names to concentrations from the NMR sample metadata.
"""
function sampleconcentrations(spec::NMRData)
    components = sample(spec, :sample, :components)
    return Dict(component["name"] => component["concentration_or_amount"]
                for component in components)
end
sampleconcentrations(expt::AbstractExperiment) = sampleconcentrations(expt.spec)

"""
    field_label(expt) -> Symbol

Convert a magnetic field strength in Tesla to a Symbol for use as a ComponentArray key.
E.g. `field_label(14.1)` → `:14p1T`.
"""
field_label(field_teslas::Float64) = Symbol(replace(string(field_teslas), "." => "p") * "T")
field_label(expt::AbstractExperiment) = field_label(expt.field_teslas)

include("expt-r1.jl")
include("expt-cest.jl")
