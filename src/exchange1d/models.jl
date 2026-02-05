# Exchange model definitions

"""
    AbstractModel

Abstract type for chemical exchange models.
"""
abstract type AbstractModel end

"""
    NoExchange <: AbstractModel

No exchange model - simple relaxation only, single state.
Used as null model for statistical comparison.
"""
struct NoExchange <: AbstractModel end

"""
    TwoState <: AbstractModel

Two-state exchange model: A ⇌ B

Exchange is characterized by:
- kex: exchange rate (kex = kAB + kBA)
- pB: population of minor state B (pA = 1 - pB)
"""
struct TwoState <: AbstractModel end

"""
    nstates(model::AbstractModel)

Return the number of states in the exchange model.
"""
nstates(::NoExchange) = 1
nstates(::TwoState) = 2

"""
    calculate_p0(model::AbstractModel, modelpars)

Calculate equilibrium populations for each state.
Returns a vector of length nstates(model).
"""
function calculate_p0(::NoExchange, modelpars)
    [1.0]
end

function calculate_p0(::TwoState, modelpars)
    pB = value(modelpars.pB)
    [1 - pB, pB]
end

"""
    calculate_K(model::AbstractModel, modelpars)

Calculate the exchange rate matrix K.
Returns an nstates × nstates matrix where K[i,j] is the rate from state j to state i.
Rows sum to zero (conservation of magnetization).
"""
function calculate_K(::NoExchange, modelpars)
    zeros(1, 1)
end

function calculate_K(::TwoState, modelpars)
    kex = value(modelpars.kex)
    pB = value(modelpars.pB)

    # Rate constants: kAB = kex * pB, kBA = kex * (1 - pB)
    # At equilibrium: pA * kAB = pB * kBA
    kAB = kex * pB
    kBA = kex * (1 - pB)

    # K matrix: dM/dt includes K*M
    # K[i,j] = rate from j to i (off-diagonal)
    # K[i,i] = -sum of rates out of i (diagonal)
    [-kAB  kBA
      kAB -kBA]
end

"""
    default_modelpars(model::AbstractModel)

Return default model parameters as a NamedTuple of Parameters.
"""
function default_modelpars(::NoExchange)
    NamedTuple()
end

function default_modelpars(::TwoState)
    (
        kex = Parameter(500.0; transform=LOG_TRANSFORM, bounds=(1.0, 1e6)),
        pB = Parameter(0.05; bounds=(0.001, 0.5)),
    )
end

"""
    modelname(model::AbstractModel)

Return a human-readable name for the model.
"""
modelname(::NoExchange) = "No Exchange"
modelname(::TwoState) = "Two-State Exchange"
