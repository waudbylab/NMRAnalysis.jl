struct InducedFitModel <: AbstractModel
    moleculemap::Dict{Symbol,String}
    concentrations::Dict{String,Float64}
end
InducedFitModel() = InducedFitModel(Dict{Symbol,String}(), Dict{String,Float64}())
InducedFitModel(moleculemap::Dict{Symbol,String}) =
    InducedFitModel(moleculemap, Dict{String,Float64}())

modelname(::InducedFitModel) = "Induced fit binding (Kd, koff, kclose, kopen)"
modelorder(::InducedFitModel) = 6
nstates(::InducedFitModel) = 3
states(::InducedFitModel) = ["free", "open", "closed"]
nmolecules(::InducedFitModel) = 2
molecules(::InducedFitModel) = Dict(:A => "observed", :X => "titrant")
function defaultparams(::InducedFitModel)
    return ComponentArray(; logKd=log(10.0), logkoff=log(1000.0),
                          logkclose=log(1000.0), logkopen=log(10.0))
end

"""
    inducedfitequilibrium(model::InducedFitModel, params, expt) -> (Afree, B, C)

Equilibrium concentrations of free observed species (`Afree`), the initial
("open") complex (`B`) and the conformationally-selected ("closed") complex
(`C`) for a two-step induced fit mechanism

    A + L ⇌ B ⇌ C

with dissociation constant `Kd` and off-rate `koff` for the bimolecular step,
and forward/backward rates `kclose`/`kopen` for the unimolecular
conformational step. Solved analytically from the ligand/observed-species
mass balance equations (Waudby, TITAN `bmInducedFit`). Shared by
`populations` and `exchangematrix`.
"""
function inducedfitequilibrium(model::InducedFitModel, params, expt)
    Kd = exp(params.model.logKd)
    kclose = exp(params.model.logkclose)
    kopen = exp(params.model.logkopen)

    A0, X0 = modelconcentrations(model, expt)

    ksum = kclose + kopen
    root = sqrt((Kd * kopen + ksum * (X0 + A0))^2 - 4 * ksum^2 * X0 * A0)

    Afree = (-Kd * kopen - ksum * (X0 - A0) + root) / (2ksum)
    B = kopen * (kclose * (X0 + A0) + kopen * (Kd + X0 + A0) - root) / (2ksum^2)
    C = kclose * (kclose * (X0 + A0) + kopen * (Kd + X0 + A0) - root) / (2ksum^2)

    return Afree, B, C
end

"""
    populations(model::InducedFitModel, params, expt) -> [pA, pB, pC]

Equilibrium populations of free (A), open-complex (B) and closed-complex (C)
species, normalised to their sum (rather than to `A0` directly) to stay exact
under floating-point roundoff.
"""
function populations(model::InducedFitModel, params, expt)
    Afree, B, C = inducedfitequilibrium(model, params, expt)
    total = Afree + B + C
    return [Afree, B, C] ./ total
end

"""
    exchangematrix(model::InducedFitModel, params, expt) -> Matrix{Float64}

Build the 3×3 kinetic exchange matrix for the induced fit mechanism
`A + L ⇌ B ⇌ C`. There is no direct A ⇌ C exchange. The pseudo-first-order
on-rate `A → B` is `kon * Lfree`, with free ligand concentration obtained
from the ligand mass balance `Lfree = X0 - (B + C)`.
"""
function exchangematrix(model::InducedFitModel, params, expt)
    Kd = exp(params.model.logKd)
    koff = exp(params.model.logkoff)
    kclose = exp(params.model.logkclose)
    kopen = exp(params.model.logkopen)
    kon = koff / Kd

    Afree, B, C = inducedfitequilibrium(model, params, expt)
    A0, X0 = modelconcentrations(model, expt)
    Lfree = X0 - (B + C)

    kAB = kon * Lfree
    kBA = koff

    return [-kAB kBA 0.0;
            kAB -kBA-kclose kopen;
            0.0 kclose -kopen]
end

"""
    modelconcentrations(model::InducedFitModel, expt) -> (A0, X0)

Look up total concentrations of the observed species and binding partner
from the experiment sample concentrations dict, using the model's `moleculemap`
to translate role symbols (:A, :X) to molecule names. Falls back to
`model.concentrations` (populated interactively) when sample metadata is
unavailable.
"""
function modelconcentrations(model::InducedFitModel, expt)
    A0 = moleculeconcentration(model, expt, :A)
    X0 = moleculeconcentration(model, expt, :X)
    return A0, X0
end
