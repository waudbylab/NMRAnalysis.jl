struct TwoStateBindingModel <: AbstractModel
    moleculemap::Dict{Symbol,String}
end
TwoStateBindingModel() = TwoStateBindingModel(Dict{Symbol,String}())

modelname(::TwoStateBindingModel) = "Two-state binding"
nstates(::TwoStateBindingModel) = 2
states(::TwoStateBindingModel) = ["free", "bound"]
nmolecules(::TwoStateBindingModel) = 2
molecules(::TwoStateBindingModel) = Dict(:A => "observed", :X => "titrant")
defaultparams(::TwoStateBindingModel) = ComponentArray(; Kd=100.0, koff=5000.0)

"""
    exchangematrix(model::TwoStateBindingModel, params, expt) -> Matrix{Float64}

Build the 2×2 kinetic exchange matrix for two-state binding.
Populations are derived from Kd and the total concentrations of the
observed species (:A) and binding partner (:X) via the quadratic
binding equation. Molecule names are looked up in `sampleconcentrations`
using `model.moleculemap`.
"""
function exchangematrix(model::TwoStateBindingModel, params, expt)
    pA, pB = populations(model, params, expt)
    koff = params.model.koff
    kon_eff = koff * pB / pA   # effective pseudo-first-order on-rate

    return [-kon_eff koff;
            kon_eff -koff]
end

function populations(model::TwoStateBindingModel, params, expt)
    Kd = params.model.Kd
    A0, X0 = modelconcentrations(model, expt)
    Xfree = 0.5 * (X0 - A0 - Kd + sqrt((Kd + A0 + X0)^2 - 4 * A0 * X0))
    Abound = X0 - Xfree
    pB = Abound / A0
    pA = 1 - pB
    return (pA, pB)
end

"""
    modelconcentrations(model::TwoStateBindingModel, expt) -> (A0, X0)

Look up total concentrations of the observed species and binding partner
from the experiment sample concentrations dict, using the model's `moleculemap`
to translate role symbols (:A, :X) to molecule names.
"""
function modelconcentrations(model::TwoStateBindingModel, expt)
    sc = sampleconcentrations(expt)
    A0 = sc[model.moleculemap[:A]]
    X0 = sc[model.moleculemap[:X]]
    return A0, X0
end
