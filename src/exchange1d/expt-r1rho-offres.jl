struct R1rhoOffResExperiment <: AbstractExperiment
    spec::Any
    field_teslas::Float64
    sampleconcentrations::Dict{String,Float64}
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
    offsets_ppm::Vector{Float64}  # offsets for each data point (ppm)
    νSL::Float64      # νSL for each data point (Hz)
    TSL::Vector{Float64}      # TSL for each data point (s)
end

function R1rhoOffResExperiment(filename)
    @info "Loading off-resonance R1rho experiment from $filename"

    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("R1rho experiment $filename must have annotations"))

    field_teslas = 2π * metadata(spec, 1, :bf) /
                   gyromagneticratio(metadata(spec, 1, :nucleus))
    field_teslas = round(field_teslas; digits=2)

    "1d" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not a 1D experiment"))
    "r1rho" in annotations(spec, :experiment_type) ||
        throw(ArgumentError("Experiment $filename is not an R1rho experiment"))
    "off_resonance" in annotations(spec, :features) ||
        throw(ArgumentError("Experiment $filename is not an off-resonance R1rho experiment"))

    # Read offsets
    offsets = annotations(spec, :r1rho, :offset)
    offsets_ppm = ppm(offsets, dims(spec, F1Dim))

    # Read spin-lock power (single value for off-resonance)
    power = annotations(spec, :r1rho, :power)
    nuc = nucleus(annotations(spec, :r1rho, :channel))
    refpulse, refpower = referencepulse(spec, nuc)
    νSL = hz(power, refpower, refpulse, 90)[1]

    # Read relaxation times
    TSL = annotations(spec, :r1rho, :duration)

    observed_intensities = zeros(length(offsets_ppm)) .± 0.0
    predicted_intensities = zeros(length(offsets_ppm))

    return R1rhoOffResExperiment(spec, field_teslas, sampleconcentrations(spec),
                                 observed_intensities, predicted_intensities,
                                 offsets_ppm, νSL, TSL)
end

function default_spin_params(expt::R1rhoOffResExperiment, nstates)
    fl = field_label(expt)
    return [:delta => fill(expt.spec[1, :offsetppm], nstates),
            Symbol("R2_", fl) => fill(10.0, nstates),
            Symbol("R1_", fl) => [1.5]]
end

function default_nuisance_params(expt::R1rhoOffResExperiment)
    return Pair{Symbol,Any}[]
end

function integrate!(expt::R1rhoOffResExperiment, peakppm, noiseppm, ppmwidth)
    spec = expt.spec

    # integrate noise region
    noiseselector = (noiseppm - ppmwidth / 2) .. (noiseppm + ppmwidth / 2)
    noise = std(vec(data(sum(spec[noiseselector, :, :]; dims=F1Dim))))

    # integrate signal region
    signalselector = (peakppm - ppmwidth / 2) .. (peakppm + ppmwidth / 2)
    integrals = sum(spec[signalselector, :, :]; dims=F1Dim)

    # normalise by max absolute value
    scale = maximum(abs, integrals)
    noise /= scale
    integrals /= scale

    # fit exponential decay for each offset
    expdecay(t, p) = @. p[1] * exp(-p[2] * t)
    p0 = [1.0, 5.0]

    for i in 1:length(expt.offsets_ppm)
        y = vec(data(integrals[1, i, :]))
        fitres = curve_fit(expdecay, expt.TSL, y, p0)
        R = coef(fitres)[2] ± stderror(fitres)[2]
        expt.observed_intensities[i] = R
    end
end

function simulate!(expt::R1rhoOffResExperiment, model, params)
    for k in 1:length(expt.offsets_ppm)
        L = liouvillian(model, params, expt,
                        expt.offsets_ppm[k],
                        expt.νSL)

        # evals = real.(eigen(L).values)
        # neg_evals = evals[evals .< 0]
        # if isempty(neg_evals)
        #     expt.predicted_intensities[k] = 1000.0
        #     continue
        # end
        # R1rho = -maximum(neg_evals)

        # Compute R1rho from inverse of trace of inverse L (Koss)
        R1rho = -1/tr(inv(L))

        expt.predicted_intensities[k] = R1rho
    end
end

function plot_result(expt::R1rhoOffResExperiment, fit_result; kwargs...)
    yobs = expt.observed_intensities
    ypred = expt.predicted_intensities
    x = expt.offsets_ppm
    sortidx = sortperm(x)

    params = fit_result.params
    params_value = fit_result.params_value

    p1 = scatter(expt.offsets_ppm, yobs;
                 xlabel="Spin-lock offset / ppm",
                 ylabel="R₁ρ / s⁻¹",
                 frame=:box, legend=nothing, grid=nothing, kwargs...)
    plot!(p1, expt.offsets_ppm[sortidx], ypred[sortidx]; lw=2)
    vline!(p1, params_value.spin.delta; ls=:dash, label="peak positions")

    wres = (Measurements.value.(yobs) .- ypred) ./ Measurements.uncertainty.(yobs)
    p2 = scatter(expt.offsets_ppm, wres;
                 xlabel="Spin-lock offset / ppm",
                 ylabel="Residual / σ",
                 frame=:box, legend=nothing, kwargs...)
    hspan!(p2, [-2, 2]; color=:limegreen, alpha=0.3, lw=0, la=0, primary=false)
    hspan!(p2, [-1, 1]; color=:limegreen, alpha=0.5, lw=0, la=0, primary=false)
    hline!(p2, [0]; color=:black, lw=0.5, grid=nothing, primary=false)
    vline!(p2, params_value.spin.delta; ls=:dash, label="peak positions", c=3)
    ylims!(p2, -maximum(abs, wres) * 1.2, maximum(abs, wres) * 1.2)

    return plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]))
end