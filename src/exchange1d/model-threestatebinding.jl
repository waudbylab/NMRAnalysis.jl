struct ThreeStateBindingModel <: AbstractModel
    moleculemap::Dict{Symbol,String}
end
ThreeStateBindingModel() = ThreeStateBindingModel(Dict{Symbol,String}())

modelname(::ThreeStateBindingModel) = "Three-state binding"
nstates(::ThreeStateBindingModel) = 3
states(::ThreeStateBindingModel) = ["free", "bound1", "bound2"]
nmolecules(::ThreeStateBindingModel) = 2
molecules(::ThreeStateBindingModel) = Dict(:A => "observed", :X => "titrant")
function defaultparams(::ThreeStateBindingModel)
    return ComponentArray(; Kd1=100.0, koff1=1000.0, Kd2=500.0, koff2=5000.0)
end

"""
    populations(model::ThreeStateBindingModel, params, expt)

Calculate equilibrium populations for three-state parallel binding:

    A + L ⇌ B  (Kd1, koff1)
    A + L ⇌ C  (Kd2, koff2)

Free ligand concentration is solved analytically from the quadratic
binding equation for two parallel binding sites (n1=n2=1).
"""

function populations(model::ThreeStateBindingModel, params, expt)
    Kd1 = params.model.Kd1
    koff1 = params.model.koff1
    Kd2 = params.model.Kd2
    koff2 = params.model.koff2

    A0, X0 = modelconcentrations(model, expt)

    # Free ligand from quadratic solution (n1=n2=1, parallel sites)
    # From TITAN bmGeneralThreeStateParallel
    Lfree = real(
        (Kd2*(X0 - A0) + Kd1*(X0 - (Kd2 + A0)) +
         sqrt(4*Kd1*Kd2*(Kd1 + Kd2)*A0 +
              (Kd1*X0 + Kd2*X0 + Kd1*Kd2 - (Kd1 + Kd2)*A0)^2)) /
        (2*(Kd1 + Kd2)))

    kon1 = koff1 / Kd1
    kon2 = koff2 / Kd2
    kab = kon1 * Lfree
    kba = koff1
    kac = kon2 * Lfree
    kca = koff2

    denom = kba*kca + kca*kab + kba*kac
    pA = kba*kca / denom
    pB = kca*kab / denom
    pC = kba*kac / denom

    return [pA, pB, pC]
end

"""
    exchangematrix(model::ThreeStateBindingModel, params, expt)

Build the 3×3 exchange matrix for three-state parallel binding.
B and C do not interconvert directly (kbc = kcb = 0).
"""
function exchangematrix(model::ThreeStateBindingModel, params, expt)
    Kd1 = params.model.Kd1
    koff1 = params.model.koff1
    Kd2 = params.model.Kd2
    koff2 = params.model.koff2

    A0, X0 = modelconcentrations(model, expt)

    Lfree = real(
        (Kd2*(X0 - A0) + Kd1*(X0 - (Kd2 + A0)) +
         sqrt(4*Kd1*Kd2*(Kd1 + Kd2)*A0 +
              (Kd1*X0 + Kd2*X0 + Kd1*Kd2 - (Kd1 + Kd2)*A0)^2)) /
        (2*(Kd1 + Kd2)))

    kon1 = koff1 / Kd1
    kon2 = koff2 / Kd2
    kab = kon1 * Lfree
    kba = koff1
    kac = kon2 * Lfree
    kca = koff2

    return [-kab-kac  kba   kca;
             kab     -kba   0.0;
             kac      0.0  -kca]
end


"""
    modelconcentrations(model::TwoStateBindingModel, expt) -> (A0, X0)

Look up total concentrations of the observed species and binding partner
from the experiment sample concentrations dict, using the model's `moleculemap`
to translate role symbols (:A, :X) to molecule names.
"""
function modelconcentrations(model::ThreeStateBindingModel, expt)
    sc = sampleconcentrations(expt)
    A0 = sc[model.moleculemap[:A]]
    X0 = sc[model.moleculemap[:X]]
    return A0, X0
end
