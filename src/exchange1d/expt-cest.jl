"""
    CESTExperiment <: AbstractExperiment

Chemical Exchange Saturation Transfer experiment. Stores saturation offsets,
B1 field strength, saturation time, and observed/predicted intensity profiles.

Loading, integration, and simulation are not yet implemented.
"""
struct CESTExperiment <: AbstractExperiment
    spec::Any                # NMRData (untyped for testability)
    field_teslas::Float64
    sampleconcentrations::Dict{String,Float64}
    δsat::Vector{Float64}                    # saturation frequencies in ppm
    ν1::Float64                           # saturation field strength in Hz
    saturation_time::Float64                    # seconds
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
end

function CESTExperiment(filename)
    detail("Loading CEST experiment from $filename")

    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("CEST experiment $filename must have annotations for analysis"))

    field_teslas = 2π * metadata(spec, 1, :bf) /
                   gyromagneticratio(metadata(spec, 1, :nucleus))
    field_teslas = round(field_teslas; digits=2)

    "1d" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not a 1D experiment"))
    "cest" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not a CEST experiment"))

    offsets = annotations(spec, :cest, :offset)
    δsat = ppm(offsets, dims(spec, F1Dim))

    satpower = annotations(spec, :cest, :power)
    nuc = nucleus(annotations(spec, :cest, :channel))
    refpulse, refpower = referencepulse(spec, nuc)
    ν1 = hz(satpower, refpower, refpulse, 90)

    saturation_time = annotations(spec, :cest, :duration)

    observed_intensities = [0.0 ± 0.0 for _ in δsat]
    predicted_intensities = zeros(Float64, length(δsat))

    return CESTExperiment(spec, field_teslas, sampleconcentrations(spec),
                          δsat, ν1, saturation_time,
                          observed_intensities, predicted_intensities)
end

"""
    default_spin_params(expt::CESTExperiment, nstates) -> Vector{Pair{Symbol,Any}}

Return spin parameter entries needed by this CEST experiment: R1, R2 for the
experiment's field, and chemical shifts (delta).
"""
function default_spin_params(expt::CESTExperiment, nstates)
    fl = field_label(expt)
    return [:delta => fill(expt.spec[1, :offsetppm], nstates),
            Symbol("R2_", fl) => fill(10.0, nstates),
            Symbol("R1_", fl) => [1.5]]
end

"""
    default_nuisance_params(expt::CESTExperiment) -> Vector{Pair{Symbol,Any}}

Return flat nuisance parameter entries for this CEST experiment: an overall
intensity scale `I0`, tagged with the experiment type and field, e.g.
`:CEST_14p1T_I0`. Shared across all CEST experiments at the same field —
see `simulate!`.
"""
function default_nuisance_params(expt::CESTExperiment)
    tag = Symbol("CEST_", field_label(expt))
    return [Symbol(tag, "_I0") => 1.0]
end

function integrate!(expt::CESTExperiment, peakppm, noiseppm, ppmwidth)
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
    simulate!(expt::CESTExperiment, model, params)

Simulate the CEST saturation profile, then scale it by the fitted `I0`
intensity (`params.nuisance.CEST_<field>_I0`) to match the observed,
integration-normalised intensities.
"""
function simulate!(expt::CESTExperiment, model, params)
    n = length(expt.δsat)
    N = nstates(model)
    T = expt.saturation_time
    p0 = populations(model, params, expt)
    M0 = zeros(3N + 1)
    for i in 1:N
        M0[3(i - 1) + 3] = p0[i]  # Mz initialised to population
    end
    M0[end] = 1.0  # augmented state for constant term

    for i in 1:n
        L = liouvillian_inhom(model, params, expt, expt.δsat[i], expt.ν1)
        M = exp(L * T) * M0
        expt.predicted_intensities[i] = sum(M[3:3:end])  # sum of Mz across states
    end

    tag = Symbol("CEST_", field_label(expt))
    I0 = params.nuisance[Symbol(tag, "_I0")]
    return expt.predicted_intensities .*= I0
end

"""
    experimentinfo(expt::CESTExperiment) -> Vector{Pair{String,String}}

Acquisition parameters worth recording alongside a saved fit (see
`_save_results`): saturation field strength, saturation time, and the
saturation offset range.
"""
function experimentinfo(expt::CESTExperiment)
    return ["Type" => "CEST",
            "Field" => _format_field(expt.field_teslas),
            "Saturation power (ν1)" => "$(round(expt.ν1; digits=1)) Hz",
            "Saturation time" => "$(round(expt.saturation_time * 1000; digits=1)) ms",
            "Saturation offsets" => "$(length(expt.δsat)) points, " *
                                    "$(round(minimum(expt.δsat); digits=3)) to " *
                                    "$(round(maximum(expt.δsat); digits=3)) ppm"]
end

function plot_result(expt::CESTExperiment, fit_result; kwargs...)
    x = expt.δsat
    sortidx = sortperm(x)  # saturation offsets may not be acquired in order

    params = fit_result.params
    params_value = fit_result.params_value

    # normalise by the fitted I0 so the baseline sits at 1.0, rather than at
    # whatever the single noisiest data point happened to integrate to
    tag = Symbol("CEST_", field_label(expt))
    I0 = params_value.nuisance[Symbol(tag, "_I0")]
    yobs = expt.observed_intensities ./ I0
    ypred = expt.predicted_intensities ./ I0

    p1 = plot(; frame=:box, legend=nothing,
              #   xlabel="Saturation frequency (ppm)",
              ylabel="Normalised intensity",
              title="CEST ($(Int(round(expt.ν1, digits=0))) Hz, $(Int(round(expt.saturation_time * 1000, digits=0))) ms)",
              grid=nothing,
              kwargs...,
              xflip=true)

    vline!(p1, params_value.spin.delta; ls=:dash, label="peak positions")
    scatter!(p1, x[sortidx], yobs[sortidx]; label="observed", ms=3)
    plot!(p1, x[sortidx], ypred[sortidx]; label="fit", lw=1.5)
    hline!(p1, [0.0]; primary=false, color=:black, lw=0.5)

    p2 = plot(; frame=:box,
              xlabel="Saturation frequency (ppm)",
              ylabel="Residual / σ",
              title="",
              grid=nothing,
              legend=nothing,
              kwargs...,
              xflip=true)
    wres = (Measurements.value.(yobs) .- ypred) ./ Measurements.uncertainty.(yobs)
    vline!(p2, params_value.spin.delta; ls=:dash, label="peak positions", c=3)
    hspan!(p2, [-2, 2]; color=:limegreen, alpha=0.3, lw=0, la=0, primary=false)
    hspan!(p2, [-1, 1]; color=:limegreen, alpha=0.5, lw=0, la=0, primary=false)
    hline!(p2, [0]; color=:black, lw=0.5, primary=false)
    scatter!(p2, x[sortidx], wres[sortidx]; ms=3)
    ylims!(p2, -maximum(abs, wres) * 1.2, maximum(abs, wres) * 1.2)

    plt = plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]), link=:x)
    return plt
end