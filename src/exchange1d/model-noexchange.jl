struct NoExchangeModel <: AbstractModel
end

modelname(::NoExchangeModel) = "No exchange"
nstates(::NoExchangeModel) = 1
nmolecules(::NoExchangeModel) = 1
molecules(::NoExchangeModel) = Dict(:X => "observed")
default_params(::NoExchangeModel) = ComponentArray()

exchange_matrix(::NoExchangeModel, params, sampleconcentrations) = zeros(1, 1)
populations(::NoExchangeModel, params, sampleconcentrations) = [1.0]
