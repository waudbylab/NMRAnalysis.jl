function componentnames(spec::NMRData)
    components = sample(spec, :sample, :components)
    isnothing(components) && return String[]
    return [component["name"] for component in components if haskey(component, "name")]
end
componentnames(expt::AbstractExperiment) = componentnames(expt.spec)

"""
    sampleconcentrations(spec::NMRData) -> Dict{String,Float64}

Extract a Dict mapping molecule names to concentrations from the NMR sample metadata.
"""
function sampleconcentrations(spec::NMRData)
    components = sample(spec, :sample, :components)
    isnothing(components) && return Dict{String,Float64}()
    result = Dict{String,Float64}()
    for component in components
        name = get(component, "name", nothing)
        conc = get(component, "concentration_or_amount", nothing)
        isnothing(name) && continue
        if isnothing(conc)
            @warn "Component \"$name\" has no concentration_or_amount defined — skipping"
            continue
        end
        result[name] = conc
    end
    return result
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
include("expt-r1rho-onres.jl")
include("expt-r1rho-offres.jl")

"""
    load_experiment(filename) -> AbstractExperiment

Load an NMR experiment file and return the appropriate concrete experiment type
based on its annotations.

Dispatches on `annotations(spec, :experiment_type)` and `annotations(spec, :features)`:
- `"relaxation"` + `"R1"` → `R1Experiment`
- `"saturation_transfer"` + `"CEST"` → `CESTExperiment` (not yet implemented)
"""
function load_experiment(filename)
    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("$filename has no annotations — cannot classify experiment"))

    types = annotations(spec, :experiment_type)
    features = annotations(spec, :features)

    if "relaxation" in types && "R1" in features
        return R1Experiment(filename)
    elseif "cest" in types
        return CESTExperiment(filename)
    elseif "r1rho" in types && "on_resonance" in features
        return R1rhoOnResExperiment(filename)
    elseif "r1rho" in types && "off_resonance" in features
        return R1rhoOffResExperiment(filename)
    else
        throw(ArgumentError("Cannot classify experiment $filename " *
                            "(types=$types, features=$features)"))
    end
end

"""
    ExchangeProblem(filenames::Vector{String}, model::AbstractModel)

Construct an ExchangeProblem by loading experiments from filenames.
Each file is classified and loaded via `load_experiment`.
"""
function ExchangeProblem(filenames::Vector{String}, model::AbstractModel)
    experiments = AbstractExperiment[load_experiment(f) for f in filenames]
    return ExchangeProblem(experiments, model)
end
