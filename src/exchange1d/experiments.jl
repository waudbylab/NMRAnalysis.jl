data(expt::AbstractExperiment) = expt.spec
concentration(expt::AbstractExperiment, component::String) = expt.concentrations[component]

function componentnames(expt::AbstractExperiment)
    components = sample(data(expt), :sample, :components)
    return [component["name"] for component in components]
end

function concentrationsdict(expt::AbstractExperiment)
    components = sample(data(expt), :sample, :components)
    return Dict(component["name"] => component["concentration"] for component in components)
end

include("expt-r1.jl")
