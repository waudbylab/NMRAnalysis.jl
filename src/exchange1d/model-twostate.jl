struct TwoStateModel <: AbstractModel
end

modelname(::TwoStateModel) = "Two-state exchange"
nstates(::TwoStateModel) = 2
states(::TwoStateModel) = ["A", "B"]
nmolecules(::TwoStateModel) = 1
molecules(::TwoStateModel) = Dict(:A => "observed")
function defaultparams(::TwoStateModel)
    pB = 0.05
    return ComponentArray(; logkex=log(1000.0), dGB=log((1 - pB) / pB))
end

"""
    exchangematrix(model::TwoStateModel, params, expt) -> Matrix{Float64}

Build the 2×2 kinetic exchange matrix for two-state intramolecular exchange.
K[i,j] is the rate from state j to state i; column sums are zero.
"""
function exchangematrix(model::TwoStateModel, params, expt)
    kex = exp(params.model.logkex)
    pA, pB = populations(model, params, expt)
    return [-kex*pB kex*pA;
            kex*pB -kex*pA]
end

"""
    populations(::TwoStateModel, params, expt) -> [pA, pB]

Population of state B is derived from the free energy difference
`dGB = G_B - G_A` (in units of RT): `pB = exp(-dGB)/(1+exp(-dGB))`, rearranged
as `1/(1+exp(dGB))` to avoid overflow for very negative `dGB`. This internal
`dGB` representation is just a fitting convenience (it keeps `pB` in (0, 1)
for any real value, unlike a raw population fraction) — the user only ever
sees/enters `pB` itself (see `_isdgparam` in interface.jl).
"""
function populations(::TwoStateModel, params, expt)
    dGB = params.model.dGB
    pB = 1 / (1 + exp(dGB))
    return [1 - pB, pB]
end
