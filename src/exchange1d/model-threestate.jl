struct ThreeStateModel <: AbstractModel
end

modelname(::ThreeStateModel) = "Three-state exchange"
nstates(::ThreeStateModel) = 3
states(::ThreeStateModel) = ["A", "B", "C"]
nmolecules(::ThreeStateModel) = 1
molecules(::ThreeStateModel) = Dict(:A => "observed")
function defaultparams(::ThreeStateModel)
    pA, pB, pC = 0.97, 0.02, 0.01
    return ComponentArray(; logkoffB=log(1000.0), dGB=log(pA / pB),
                          logkoffC=log(500.0), dGC=log(pA / pC))
end

function exchangematrix(model::ThreeStateModel, params, expt)
    pA, pB, pC = populations(model, params, expt)
    koffB = exp(params.model.logkoffB)
    koffC = exp(params.model.logkoffC)
    konB = koffB * pB / pA
    konC = koffC * pC / pA
    return [-konB-konC koffB koffC;
            konB -koffB 0.0;
            konC 0.0 -koffC]
end

"""
    populations(::ThreeStateModel, params, expt) -> [pA, pB, pC]

Populations are derived from the free energy differences `dGB = G_B - G_A`
and `dGC = G_C - G_A` (in units of RT), relative to reference state A. This
keeps all three populations in (0, 1) summing to 1 for any real `dGB`,
`dGC`. This internal `dG` representation is just a fitting convenience —
the user only ever sees/enters `pB`/`pC` themselves (see `_isdgparam` in
interface.jl).
"""
function populations(::ThreeStateModel, params, expt)
    dGB = params.model.dGB
    dGC = params.model.dGC
    denom = 1 + exp(-dGB) + exp(-dGC)
    pA = 1 / denom
    pB = exp(-dGB) / denom
    pC = exp(-dGC) / denom
    return [pA, pB, pC]
end
