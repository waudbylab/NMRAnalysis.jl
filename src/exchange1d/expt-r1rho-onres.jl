struct R1rhoOnResExperiment <: AbstractExperiment
    spec::Any
    field_teslas::Float64
    sampleconcentrations::Dict{String,Float64}
    observed_intensities::Vector{Measurement{Float64}}  # observed (fitted) relaxation rates vs vSL
    predicted_intensities::Vector{Float64}              # predicted relaxation rates vs vSL

    νSL::Vector{Float64}      # νSL for each data point (Hz)
    TSL::Vector{Float64}      # TSL for each data point (s)
end

function R1rhoOnResExperiment(filename)
    detail("Loading on-resonance R1rho experiment from $filename")

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
    "on_resonance" in annotations(spec, :features) ||
        throw(ArgumentError("Experiment $filename is not an on-resonance R1rho experiment"))

    # Read spin-lock powers and convert to Hz
    powers = annotations(spec, :r1rho, :power)
    nuc = nucleus(annotations(spec, :r1rho, :channel))
    refpulse, refpower = referencepulse(spec, nuc)
    νSL = hz.(powers, refpower, refpulse, 90)

    # Read relaxation times
    TSL = annotations(spec, :r1rho, :duration)

    observed_intensities = zeros(length(νSL)) .± 0.0
    predicted_intensities = zeros(length(νSL))

    return R1rhoOnResExperiment(spec, field_teslas, sampleconcentrations(spec),
                                observed_intensities, predicted_intensities, νSL, TSL)
end

function default_spin_params(expt::R1rhoOnResExperiment, nstates)
    fl = field_label(expt)
    δs = fill(expt.spec[1, :offsetppm], nstates)
    δs[2:end] .+= randn(nstates - 1)
    return [:delta => δs,
            Symbol("R2_", fl) => fill(10.0, nstates),
            Symbol("R1_", fl) => [1.5]]
end

function default_nuisance_params(expt::R1rhoOnResExperiment)
    return Pair{Symbol,Any}[]
end

# integrate peak and fit to exponential to determine R1rho rates
function integrate!(expt::R1rhoOnResExperiment, peakppm, noiseppm, ppmwidth)
    spec = expt.spec

    # integrate noise region
    noiseselector = (noiseppm - ppmwidth / 2) .. (noiseppm + ppmwidth / 2)
    n = sum(spec[noiseselector, :, :]; dims=F1Dim)

    noise = vec(data(n))
    noise = std(noise)

    # integrate signal region
    signalselector = (peakppm - ppmwidth / 2) .. (peakppm + ppmwidth / 2)
    integrals = sum(spec[signalselector, :, :]; dims=F1Dim)

    # normalise by max absolute value
    scale = maximum(abs, integrals)
    noise /= scale
    integrals /= scale

    # fit to exponential decays
    expdecay(t, p) = @. p[1] * exp(-p[2] * t)
    p0 = [1.0, 20.0]

    for i in 1:length(expt.νSL)
        y = vec(data(integrals[1, i, :]))
        fitres = curve_fit(expdecay, expt.TSL, y, p0)
        R = coef(fitres)[2] ± stderror(fitres)[2]
        expt.observed_intensities[i] = R
    end
end

function simulate!(expt::R1rhoOnResExperiment, model, params)
    fl = field_label(expt)
    spinlock_ppm = params.spin.delta[1]

    for k in 1:length(expt.νSL)
        L = liouvillian(model, params, expt,
                        spinlock_ppm,
                        expt.νSL[k])
        # Exact R1rho from least negative eigenvalue of L
        # evals = real.(eigen(L).values)
        # R1rho = -maximum(evals[evals .< 0])

        # Compute R1rho from inverse of trace of inverse L (Koss)
        R1rho = -1/tr(inv(L))

        expt.predicted_intensities[k] = R1rho
    end
end

"""
    experimentinfo(expt::R1rhoOnResExperiment) -> Vector{Pair{String,String}}

Acquisition parameters worth recording alongside a saved fit (see
`_save_results`): field strength and spin-lock strength range.
"""
function experimentinfo(expt::R1rhoOnResExperiment)
    return ["Type" => "On-resonance R1ρ",
            "Field" => _format_field(expt.field_teslas),
            "Spin-lock strengths" =>
                "$(length(expt.νSL)) points, " *
                "$(round(minimum(expt.νSL); digits=1)) to " *
                "$(round(maximum(expt.νSL); digits=1)) Hz"]
end

function plotresult!(gl, expt::R1rhoOnResExperiment, fitresult; axiskw=(;))
    yobs = expt.observed_intensities
    ypred = expt.predicted_intensities
    x = expt.νSL
    sortidx = sortperm(x)
    wres = (Measurements.value.(yobs) .- ypred) ./ Measurements.uncertainty.(yobs)

    ax1, ax2 = resultpanels!(gl; xlabel="Spin-lock strength / Hz", ylabel="R₁ρ / s⁻¹",
                             title="On-resonance R₁ρ", axiskw=axiskw)

    hlines!(ax1, [0.0]; color=:black, linewidth=0.5)
    measured!(ax1, x, yobs; color=palettecolor(1))
    lines!(ax1, x[sortidx], ypred[sortidx]; color=palettecolor(2), linewidth=1.5)

    residualbands!(ax2)
    scatter!(ax2, x, wres; color=palettecolor(1), markersize=6)
    ylims!(ax2, -maximum(abs, wres) * 1.2, maximum(abs, wres) * 1.2)
    return ax1, ax2
end
