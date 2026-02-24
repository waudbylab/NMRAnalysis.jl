struct TwoStateBindingModel <: AbstractModel
    moleculemap::Dict{Symbol,String}
end
TwoStateBindingModel() = TwoStateBindingModel(Dict{Symbol,String}())

modelname(::TwoStateBindingModel) = "Two-state binding"
nstates(::TwoStateBindingModel) = 2
nmolecules(::TwoStateBindingModel) = 2
molecules(::TwoStateBindingModel) = Dict(:X => "observed", :Y => "binding partner")
default_params(::TwoStateBindingModel) = ComponentArray(; Kd=100.0, koff=5000.0)

"""
    exchange_matrix(model::TwoStateBindingModel, params, sampleconcentrations) -> Matrix{Float64}

Build the 2×2 kinetic exchange matrix for two-state binding.
Populations are derived from Kd and the total concentrations of the
observed species (:X) and binding partner (:Y) via the quadratic
binding equation. Molecule names are looked up in `sampleconcentrations`
using `model.moleculemap`.
"""
function exchange_matrix(model::TwoStateBindingModel, params, expt)
    Kd = params.model.Kd
    koff = params.model.koff

    Xt, Yt = _lookup_concentrations(model, sampleconcentrations(expt))
    pB = _binding_fraction(Kd, Xt, Yt)
    pA = 1 - pB
    kon_eff = koff * pB / pA   # effective pseudo-first-order on-rate

    return [-kon_eff koff;
            kon_eff -koff]
end

function populations(model::TwoStateBindingModel, params, expt)
    Kd = params.model.Kd
    Xt, Yt = _lookup_concentrations(model, sampleconcentrations(expt))
    pB = _binding_fraction(Kd, Xt, Yt)
    return [1 - pB, pB]
end

"""
    _lookup_concentrations(model::TwoStateBindingModel, sampleconcentrations) -> (Xt, Yt)

Look up total concentrations of the observed species and binding partner
from the sample concentrations dict, using the model's `moleculemap` to
translate role symbols (:X, :Y) to molecule names.
"""
function _lookup_concentrations(model::TwoStateBindingModel, sampleconcentrations)
    if isempty(model.moleculemap)
        throw(ArgumentError("TwoStateBindingModel.moleculemap must be set before use"))
    end
    Xt = sampleconcentrations[model.moleculemap[:X]]
    Yt = sampleconcentrations[model.moleculemap[:Y]]
    return Xt, Yt
end

"""
    _binding_fraction(Kd, Xt, Yt) -> Float64

Calculate the fraction bound (pB) from the dissociation constant Kd and total
concentrations Xt (observed species) and Yt (binding partner), using the
standard quadratic binding equation.
"""
function _binding_fraction(Kd, Xt, Yt)
    # Quadratic: pB² · Xt - pB · (Xt + Yt + Kd) + Yt = 0
    a = Xt
    b = -(Xt + Yt + Kd)
    c = Yt
    discriminant = b^2 - 4 * a * c
    discriminant < 0 && throw(DomainError(discriminant,
                                          "Negative discriminant in binding equation — check Kd and concentrations"))
    return (-b - sqrt(discriminant)) / (2 * a)
end
