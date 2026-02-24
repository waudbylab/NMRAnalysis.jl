struct TwoStateModel <: AbstractModel
end

modelname(::TwoStateModel) = "Two-state exchange"
nstates(::TwoStateModel) = 2
nmolecules(::TwoStateModel) = 1
molecules(::TwoStateModel) = Dict(:X => "observed")
default_params(::TwoStateModel) = ComponentArray(; kex=1000.0, pB=0.05)