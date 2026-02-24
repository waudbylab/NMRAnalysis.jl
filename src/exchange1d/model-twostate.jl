struct TwoStateModel <: AbstractModel
end

modelname(::TwoStateModel) = "Two-state exchange"
nstates(::TwoStateModel) = 2
nmolecules(::TwoStateModel) = 1
molecules(::TwoStateModel) = Dict(:X => "observed")
default_params(::TwoStateModel) = ComponentArray(; kex=1000.0, pB=0.05)

"""
    exchange_matrix(::TwoStateModel, params, sampleconcentrations) -> Matrix{Float64}

Build the 2×2 kinetic exchange matrix for two-state intramolecular exchange.
K[i,j] is the rate from state j to state i; column sums are zero.
"""
function exchange_matrix(::TwoStateModel, params, sampleconcentrations)
    kex = params.model.kex
    pB = params.model.pB
    pA = 1 - pB
    return [-kex * pB   kex * pA;
             kex * pB  -kex * pA]
end

function populations(::TwoStateModel, params, sampleconcentrations)
    pB = params.model.pB
    return [1 - pB, pB]
end
