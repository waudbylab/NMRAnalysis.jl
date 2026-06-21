struct CustomModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    initial_values::Dict{String,Float64}
end

struct ExponentialModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

struct RecoveryModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

function estimate_parameters(x::AbstractVector, y::AbstractVector, ::RecoveryModel)
    A = maximum(y)
    C = 1.0
    R = 3.0 / maximum(x)
    Dict("A" => A, "C" => C, "R" => R)
end

function estimate_parameters(x::AbstractVector, y::AbstractVector, ::ExponentialModel)
    A = maximum(abs.(y))
    R = 3.0 / maximum(x)  # Different heuristic for T1
    Dict("A" => A, "R" => R)
end

function estimate_parameters(x::AbstractVector, y::AbstractVector, model::CustomModel)
    model.initial_values
end

function ExponentialModel()
    ExponentialModel(
        (x, p) -> (@. p[1] * exp(-p[2] * x)),
        ["A", "R"],
        "Time / s"
    )
end

function RecoveryModel()
    RecoveryModel(
        (x, p) -> (@. p[1] * (1 - p[2] * exp(-p[3] * x))),
        ["A", "C", "R"],
        "Time / s"
    )
end

function CustomModel(modelfunction::String, params::Vector{Pair{String,Float64}}, xlabel="x")
    param_names = first.(params)
    CustomModel(
        modeltofunction(modelfunction, param_names),
        param_names,
        xlabel,
        Dict(params)
    )
end

function modeltofunction(expr, param_names)
    function replace_walker(ex)
        if ex isa Symbol && string(ex) in param_names
            # Found a parameter name as a symbol - replace with indexed access
            idx = findfirst(==(string(ex)), param_names)
            return :(p[$idx])
        elseif ex isa Expr
            # Recursively walk through expression
            return Expr(ex.head, map(replace_walker, ex.args)...)
        else
            # Leave other elements (numbers, etc) unchanged
            return ex
        end
    end
    parsed_expr = Meta.parse(expr)
    expr_with_params = replace_walker(parsed_expr)
    func_expr = :($(Expr(:tuple, :x, :p)) -> @. $expr_with_params)
    eval(func_expr)
end

function postfit!(peak::Peak, expt::IntensityExperiment, ::NoFitting)
    peak.postfitted[] = true
end

# Generic parametric model fitting
function postfit!(peak::Peak, expt::IntensityExperiment, model::ParametricModel)
    @debug "Post-fitting model"
    x = expt.x
    y = peak.parameters[:amp].value[]

    # Exclude skipped planes from the fit
    skip = Set(expt.skipplanes)
    keep = [i for i in eachindex(x) if i ∉ skip]
    x_fit = x[keep]
    y_fit = y[keep]

    p0 = collect(values(estimate_parameters(x_fit, y_fit, model)))
    fit = curve_fit(model.func, x_fit, y_fit, p0)
    pfit = coef(fit)
    perr = stderror(fit)

    for (i, name) in enumerate(model.param_names)
        param = peak.postparameters[Symbol(name)]
        param.value[] .= pfit[i]
        param.uncertainty[] .= perr[i]
    end

    @debug "Fitted parameters: $(peak.postparameters)"
    peak.postfitted[] = true
end

# Helper: plane indices to skip (empty for experiments that don't support skipplanes)
_skipset(::Experiment) = Set{Int}()
_skipset(expt::IntensityExperiment) = Set{Int}(expt.skipplanes)

_empty_errorbars() = Tuple{Float64,Float64,Float64}[]

# Default - just return amplitudes, separating active and skipped points
function get_model_data(peak, expt::Experiment, ::NoFitting)
    isnothing(peak) &&
        return (Point2f[], _empty_errorbars(), Point2f[], Point2f[], _empty_errorbars())

    x = expt.x
    y = peak.parameters[:amp].value[]
    err = peak.parameters[:amp].uncertainty[]
    skip = _skipset(expt)

    active  = [i for i in eachindex(x) if i ∉ skip]
    skipped = [i for i in eachindex(x) if i ∈ skip]

    obs_points    = Point2f.(x[active], y[active])
    obs_errors    = [(x[i], y[i], err[i]) for i in active]
    skip_points   = Point2f.(x[skipped], y[skipped])
    skip_errors   = [(x[i], y[i], err[i]) for i in skipped]

    return (obs_points, obs_errors, Point2f[], skip_points, skip_errors)
end

# Generic parametric model visualization
function get_model_data(peak, expt::Experiment, model::ParametricModel)
    isnothing(peak) &&
        return (Point2f[], _empty_errorbars(), Point2f[], Point2f[], _empty_errorbars())

    x = expt.x
    y = peak.parameters[:amp].value[]
    err = peak.parameters[:amp].uncertainty[]
    skip = _skipset(expt)

    active  = [i for i in eachindex(x) if i ∉ skip]
    skipped = [i for i in eachindex(x) if i ∈ skip]

    obs_points  = Point2f.(x[active], y[active])
    obs_errors  = [(x[i], y[i], err[i]) for i in active]
    skip_points = Point2f.(x[skipped], y[skipped])
    skip_errors = [(x[i], y[i], err[i]) for i in skipped]

    if peak.postfitted[]
        xpred = range(min(0.0, minimum(x)), 1.1 * maximum(x), 100)
        p = [peak.postparameters[Symbol(name)].value[][1] for name in model.param_names]
        ypred = model.func(xpred, p)
        fit_points = Point2f.(xpred, ypred)
    else
        fit_points = Point2f[]
    end

    return (obs_points, obs_errors, fit_points, skip_points, skip_errors)
end

function model_parameter_text(peak::Peak, ::NoFitting)
    ["Amplitude: $(peak.parameters[:amp].value[][1] ± peak.parameters[:amp].uncertainty[][1])"]
end

function model_parameter_text(peak::Peak, model::ParametricModel)
    map(model.param_names) do name
        param = peak.postparameters[Symbol(name)]
        "$name: $(param.value[][1] ± param.uncertainty[][1])"
    end
end

function model_info_text(::NoFitting, x::AbstractVector)
    ["Number of spectra: $(length(x))"]
end

function model_info_text(model::ParametricModel, x::AbstractVector)
    [
        "Number of points: $(length(x))",
        "$(model.xlabel) range: $(minimum(x)) - $(maximum(x))"
    ]
end