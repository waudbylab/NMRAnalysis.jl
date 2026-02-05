# Parameter type for model fitting with transforms and bounds

"""
    Parameter

Represents a fittable parameter with value, bounds, optional transform, and fixed state.

# Fields
- `initial::Float64`: Initial value for fitting
- `fitted::Union{Nothing, Float64}`: Fitted value (nothing if not yet fitted)
- `transform::Tuple{Function, Function}`: (to_fit_space, from_fit_space) transforms
- `bounds::Tuple{Float64, Float64}`: (lower, upper) bounds in original space
- `fixed::Bool`: Whether parameter is fixed during fitting
"""
struct Parameter
    initial::Float64
    fitted::Union{Nothing, Float64}
    transform::Tuple{Function, Function}
    bounds::Tuple{Float64, Float64}
    fixed::Bool
end

"""
    Parameter(initial; transform=(identity, identity), bounds=(-Inf, Inf), fixed=false)

Create a new Parameter with given initial value and optional constraints.
"""
function Parameter(initial::Real; transform=(identity, identity), bounds=(-Inf, Inf), fixed=false)
    Parameter(Float64(initial), nothing, transform, Float64.(bounds), fixed)
end

"""
    fixed(val)

Create a fixed Parameter that will not be varied during fitting.
"""
function fixed(val::Real)
    v = Float64(val)
    Parameter(v, v, (identity, identity), (v, v), true)
end

"""
    value(p::Parameter)

Get the current value of a parameter (fitted value if available, otherwise initial).
"""
value(p::Parameter) = isnothing(p.fitted) ? p.initial : p.fitted

"""
    to_fit(p::Parameter)

Transform parameter value to fitting space.
"""
to_fit(p::Parameter) = p.transform[1](value(p))

"""
    from_fit(p::Parameter, x)

Transform value from fitting space back to original space.
"""
from_fit(p::Parameter, x) = p.transform[2](x)

"""
    with_fitted(p::Parameter, fitted_value)

Return a new Parameter with the fitted value set.
"""
function with_fitted(p::Parameter, fitted_value::Real)
    Parameter(p.initial, Float64(fitted_value), p.transform, p.bounds, p.fixed)
end

# Common transforms
"""Log transform for positive parameters spanning orders of magnitude (e.g., kex, amplitude)"""
const LOG_TRANSFORM = (log, exp)
