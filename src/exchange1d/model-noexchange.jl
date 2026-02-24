struct NoExchangeModel <: AbstractModel
end

modelname(::NoExchangeModel) = "No exchange"
nstates(::NoExchangeModel) = 1
nmolecules(::NoExchangeModel) = 1
molecules(::NoExchangeModel) = Dict(:X => "observed")
default_params(::NoExchangeModel) = ComponentArray()