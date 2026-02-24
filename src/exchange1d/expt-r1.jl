struct R1Experiment <: AbstractExperiment
    spec::NMRData
    field_teslas::Float64
    concentrations::Dict{String,Float64}
    delays::Vector{Float64}
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
end

function R1Experiment(filename)
    @info "Loading R1 experiment from $filename"

    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("R1 experiment $filename must have annotations for analysis"))

    field_teslas = 2π * metadata(spec, 1, :bf) /
                   gyromagneticratio(metadata(spec, 1, :nucleus))
    field_teslas = round(field_teslas; digits=2)

    type = annotations(spec, :relaxation, :type)
    type == "R1" ||
        throw(ArgumentError("Experiment $filename is not an R1 relaxation experiment"))

    tau = annotations(spec, :relaxation, :duration)

    observed_intensities = [0.0 ± 0.0 for _ in tau]
    predicted_intensities = zeros(Float64, length(tau))

    return R1Experiment(spec, field_teslas, concentrationsdict(spec), tau,
                        observed_intensities, predicted_intensities)
end

function integrate!(expt::R1Experiment, peakppm::Float64, noiseppm::Float64,
                    ppmwidth::Float64)
    # Placeholder for integration logic
    @info "Integrating R1 experiment at peak position $peakppm with noise position $noiseppm and width $ppmwidth"

    spec = expt.spec

    # integrate regions
    noiseselector = (noiseppm - ppmwidth / 2) .. (noiseppm + ppmwidth / 2)
    noise = vec(data(sum(spec[noiseselector, :]; dims=F1Dim)))
    noise = std(noise)

    signalselector = (peakppm - ppmwidth / 2) .. (peakppm + ppmwidth / 2)
    integrals = vec(data(sum(spec[signalselector, :]; dims=F1Dim)))

    # normalise by max value
    noise /= maximum(integrals)
    integrals /= maximum(integrals)

    # update vector containing observed intensities
    return expt.observed_intensities .= integrals .± noise
end

# ir = annotations(spec, :relaxation, :model) == "inversion_recovery"
# nuc = annotations(spec, :relaxation, :channel)

# model(t, p) = ir ? p[1] * (1 .- p[3] * exp.(-t * p[2])) : p[1] * exp.(-t * p[2])
# p0 = ir ? [1.0, 2 / maximum(tau), 2.0] : [1.0, 2 / maximum(tau)]