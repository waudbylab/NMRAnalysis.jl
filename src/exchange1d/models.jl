
struct NoExchangeModel <: AbstractModel
end
modelname(::NoExchangeModel) = "No exchange"
nstates(::NoExchangeModel) = 1
nmolecules(::NoExchangeModel) = 1
molecules(::NoExchangeModel) = Dict(:X => "observed")
default_params(::NoExchangeModel) = ComponentArray()

struct TwoStateModel <: AbstractModel
end
modelname(::TwoStateModel) = "Two-state exchange"
nstates(::TwoStateModel) = 2
nmolecules(::TwoStateModel) = 1
molecules(::TwoStateModel) = Dict(:X => "observed")
default_params(::TwoStateModel) = ComponentArray(; kex=1000.0, pB=0.05)

struct TwoStateBindingModel <: AbstractModel
    moleculemap = Dict{Symbol,String}
end
modelname(::TwoStateBindingModel) = "Two-state binding"
nstates(::TwoStateBindingModel) = 2
nmolecules(::TwoStateBindingModel) = 2
molecules(::TwoStateBindingModel) = Dict(:X => "observed", :Y => "binding partner")
default_params(::TwoStateBindingModel) = ComponentArray(; Kd=100.0, koff=5000.0)