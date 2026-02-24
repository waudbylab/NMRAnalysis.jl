struct TwoStateBindingModel <: AbstractModel
    moleculemap = Dict{Symbol,String}
end

modelname(::TwoStateBindingModel) = "Two-state binding"
nstates(::TwoStateBindingModel) = 2
nmolecules(::TwoStateBindingModel) = 2
molecules(::TwoStateBindingModel) = Dict(:X => "observed", :Y => "binding partner")
default_params(::TwoStateBindingModel) = ComponentArray(; Kd=100.0, koff=5000.0)