struct R1Experiment <: AbstractExperiment
    spec::Any                # NMRData (untyped for testability)
    field_teslas::Float64
    sampleconcentrations::Dict{String,Float64}
    delays::Vector{Float64}
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
    fitting_model::Symbol   # :exponential_decay or :inversion_recovery
end

function R1Experiment(filename)
    detail("Loading R1 experiment from $filename")

    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("R1 experiment $filename must have annotations for analysis"))

    field_teslas = 2π * metadata(spec, 1, :bf) /
                   gyromagneticratio(metadata(spec, 1, :nucleus))
    field_teslas = round(field_teslas; digits=2)

    "1d" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not a 1D experiment"))
    "relaxation" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not a relaxation experiment"))
    type = annotations(spec, :relaxation, :type)
    type == "R1" ||
        throw(ArgumentError("Experiment $filename is not an R1 relaxation experiment"))

    tau = annotations(spec, :relaxation, :duration)

    fitting_model = if annotations(spec, :relaxation, :model) == "inversion_recovery"
        :inversion_recovery
    else
        :exponential_decay
    end

    observed_intensities = [0.0 ± 0.0 for _ in tau]
    predicted_intensities = zeros(Float64, length(tau))

    return R1Experiment(spec, field_teslas, sampleconcentrations(spec), tau,
                        observed_intensities, predicted_intensities, fitting_model)
end

function integrate!(expt::R1Experiment, peakppm, noiseppm, ppmwidth)
    spec = expt.spec

    # integrate noise region
    noiseselector = (noiseppm - ppmwidth / 2) .. (noiseppm + ppmwidth / 2)
    n = sum(spec[noiseselector, :]; dims=F1Dim)

    noise = vec(data(n))
    noise = std(noise)

    # integrate signal region
    signalselector = (peakppm - ppmwidth / 2) .. (peakppm + ppmwidth / 2)
    integrals = vec(data(sum(spec[signalselector, :]; dims=F1Dim)))

    # normalise by max absolute value
    scale = maximum(abs, integrals)
    noise /= scale
    integrals /= scale

    # update vector containing observed intensities
    return expt.observed_intensities .= integrals .± noise
end

"""
    default_spin_params(expt::R1Experiment, nstates) -> Vector{Pair{Symbol,Any}}

Return spin parameter entries needed by this R1 experiment: R1 (shared, length=1)
"""
function default_spin_params(expt::R1Experiment, nstates)
    fl = field_label(expt)
    return [Symbol("R1_", fl) => [1.5]]
end

"""
    default_nuisance_params(expt::R1Experiment) -> Vector{Pair{Symbol,Any}}

Return flat nuisance parameter entries for this R1 experiment, tagged with
the experiment type and field, e.g. `:R1_14p1T_I0`, `:R1_14p1T_inv_factor`.
"""
function default_nuisance_params(expt::R1Experiment)
    fl = field_label(expt)
    pairs = Pair{Symbol,Any}[Symbol("R1_", fl, "_I0") => 1.0]
    if expt.fitting_model == :inversion_recovery
        push!(pairs, Symbol("R1_", fl, "_inv_factor") => 2.0)
    end
    return pairs
end

"""
    simulate!(expt::R1Experiment, model::AbstractModel, params::ComponentArray)

Simulate predicted intensities for an R1 experiment. The `model` argument is
accepted for interface consistency but is not used (R1 decay is independent of
chemical exchange).

Reads R1 from `params.spin.R1_<field>` and nuisance parameters
(e.g. `params.nuisance.R1_14p1T_I0`) directly.
"""
function simulate!(expt::R1Experiment, ::AbstractModel, params::ComponentArray)
    fl = field_label(expt)
    R1 = params.spin[Symbol("R1_", fl)][1]

    tag = Symbol("R1_", fl)
    I0 = params.nuisance[Symbol(tag, "_I0")]

    if expt.fitting_model == :inversion_recovery
        inv_factor = params.nuisance[Symbol(tag, "_inv_factor")]
        expt.predicted_intensities .= I0 .* (1 .- inv_factor .* exp.(-expt.delays .* R1))
    else
        expt.predicted_intensities .= I0 .* exp.(-expt.delays .* R1)
    end

    return nothing
end

"""
    experimentinfo(expt::R1Experiment) -> Vector{Pair{String,String}}

Acquisition parameters worth recording alongside a saved fit (see
`_save_results`): field strength, fitting model, and relaxation delay range.
"""
function experimentinfo(expt::R1Experiment)
    return ["Type" => "R1 relaxation",
            "Field" => _format_field(expt.field_teslas),
            "Fitting model" => string(expt.fitting_model),
            "Delays" =>
                "$(length(expt.delays)) points, " *
                "$(round(minimum(expt.delays); digits=3)) to " *
                "$(round(maximum(expt.delays); digits=3)) s"]
end

"""
    plotresult!(gl, expt::R1Experiment, fitresult; axiskw=(;))

Draw an R1 experiment into `gl`: intensities against delay with the fitted curve above,
weighted residuals below.
"""
function plotresult!(gl, expt::R1Experiment, fitresult; axiskw=(;))
    param_values = fitresult.params_value
    tau = expt.delays
    yobs = expt.observed_intensities
    ypred = expt.predicted_intensities

    fl = field_label(expt)
    R1_value = param_values.spin[Symbol("R1_", fl)][1]

    # smooth fitted curve
    x = LinRange(0, maximum(tau) * 1.1, 100)
    tag = Symbol("R1_", fl)
    I0 = param_values.nuisance[Symbol(tag, "_I0")]
    yfit = if expt.fitting_model == :inversion_recovery
        inv_factor = param_values.nuisance[Symbol(tag, "_inv_factor")]
        I0 .* (1 .- inv_factor .* exp.(-x .* R1_value))
    else
        I0 .* exp.(-x .* R1_value)
    end

    wres = (Measurements.value.(yobs) .- ypred) ./ Measurements.uncertainty.(yobs)
    wrange = maximum(abs, wres) * 1.2

    ax1, ax2 = resultpanels!(gl; xlabel="Relaxation time / s", ylabel="Intensity",
                             title="R1 relaxation", axiskw=axiskw)

    hlines!(ax1, [0]; color=:black, linewidth=0.5)
    measured!(ax1, tau, yobs; color=palettecolor(1))
    lines!(ax1, x, yfit; color=palettecolor(2), linewidth=1.5)

    residualbands!(ax2)
    scatter!(ax2, tau, wres; color=palettecolor(1), markersize=6)
    ylims!(ax2, -wrange, wrange)
    return ax1, ax2
end
