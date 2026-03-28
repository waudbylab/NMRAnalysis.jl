struct ThreeStateModel <: AbstractModel
end

modelname(::ThreeStateModel) = "Three-state exchange"
nstates(::ThreeStateModel) = 3
states(::ThreeStateModel) = ["A", "B", "C"]
nmolecules(::ThreeStateModel) = 1
molecules(::ThreeStateModel) = Dict(:A => "observed")
function defaultparams(::ThreeStateModel)
    return ComponentArray(; koffB=1000.0, pB=0.02, koffC=500.0, pC=0.01)
end

function exchangematrix(::ThreeStateModel, params, expt)
    pB = params.model.pB
    pC = params.model.pC
    pA = 1 - pB - pC
    koffB = params.model.koffB
    koffC = params.model.koffC
    konB = koffB * pB / pA
    konC = koffC * pC / pA
    return [-konB-konC konB konC;
            koffB -koffB 0.0;
            koffC 0.0 -koffC]
end

function populations(::ThreeStateModel, params, expt)
    pB = params.model.pB
    pC = params.model.pC
    return [1 - pB - pC, pB, pC]
end
