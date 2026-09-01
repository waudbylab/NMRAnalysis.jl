struct NoExchangeModel <: AbstractModel
end

modelname(::NoExchangeModel) = "No exchange"
nstates(::NoExchangeModel) = 1
states(::NoExchangeModel) = ["A"]
nmolecules(::NoExchangeModel) = 1
molecules(::NoExchangeModel) = Dict(:X => "observed")
defaultparams(::NoExchangeModel) = ComponentArray()

exchangematrix(::NoExchangeModel, params, expt) = zeros(1, 1)
populations(::NoExchangeModel, params, expt) = [1.0]
