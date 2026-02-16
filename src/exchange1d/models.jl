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
    return [1.0]
end

function calculate_p0(::TwoState, modelpars)
    pB = value(modelpars.pB)
    return [1 - pB, pB]
end

"""
    calculate_K(model::AbstractModel, modelpars)

Calculate the exchange rate matrix K.
Returns an nstates × nstates matrix where K[i,j] is the rate from state j to state i.
Rows sum to zero (conservation of magnetization).
"""
function calculate_K(::NoExchange, modelpars)
    return zeros(1, 1)
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
    return [-kAB kBA
            kAB -kBA]
end

"""
    default_modelpars(model::AbstractModel)

Return default model parameters as a NamedTuple of Parameters.
"""
function default_modelpars(::NoExchange)
    return NamedTuple()
end

function default_modelpars(::TwoState)
    return (kex=Parameter(500.0; transform=LOG_TRANSFORM, bounds=(1.0, 1e6)),
            pB=Parameter(0.05; bounds=(0.001, 0.5)))
end

"""
    modelname(model::AbstractModel)

Return a human-readable name for the model.
"""
modelname(::NoExchange) = "No Exchange"
modelname(::TwoState) = "Two-State Exchange"

"""
    select_models_menu()

Display multi-select menu for model selection.
"""
function select_models_menu()
    options = ["No exchange (null model)",
               "Two-state exchange"]
    menu = MultiSelectMenu(options)
    choices = request("Select models to fit:", menu)

    models = AbstractModel[]
    1 in choices && push!(models, NoExchange())
    2 in choices && push!(models, TwoState())

    if isempty(models)
        error("No models selected")
    end

    return models
end

"""
    prompt_model_parameters(model::NoExchange)

No parameters needed for NoExchange model.
"""
function prompt_model_parameters(::NoExchange)
    return NamedTuple()
end

"""
    prompt_model_parameters(model::TwoState)

Prompt user for two-state exchange parameters (kex, pB).
"""
function prompt_model_parameters(::TwoState)
    println("Enter initial parameter estimates for two-state exchange:")

    print("  Exchange rate kex (s⁻¹) [default 500]: ")
    kex_input = readline()
    kex = isempty(strip(kex_input)) ? 500.0 : parse(Float64, kex_input)

    print("  Minor state population pB [default 0.05]: ")
    pB_input = readline()
    pB = isempty(strip(pB_input)) ? 0.05 : parse(Float64, pB_input)

    return (kex=Parameter(kex; transform=LOG_TRANSFORM, bounds=(1.0, 1e6)),
            pB=Parameter(pB; bounds=(0.001, 0.5)))
end