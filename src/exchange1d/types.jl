abstract type AbstractModel end

abstract type AbstractExperiment end

struct ExchangeProblem
    experiments::Vector{AbstractExperiment}
    model::AbstractModel
end