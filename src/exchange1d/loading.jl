"""
    LoadedData

Container for loaded experiment data before integration.
"""
struct LoadedData
    Ω::Vector{Float64}      # Offsets (rad/s)
    ω1::Vector{Float64}     # Field strengths (rad/s)
    t::Vector{Float64}      # Times (s)
    spectra::Vector{NMRData}  # 1D spectra
    bf::Float64             # Spectrometer frequency (Hz)
end

"""
    load_r1rho_data(filenames::Vector{String}; minω1=0.0, maxω1=Inf)

Load R1ρ experiments from NMR data files.

# Arguments
- `filenames`: Paths to Bruker experiment directories
- `minω1`: Minimum spin-lock field strength to include (rad/s)
- `maxω1`: Maximum spin-lock field strength to include (rad/s)

# Returns
- `LoadedData` containing conditions and spectra
"""
function load_r1rho_data(filenames::Vector{String}; minω1=0.0, maxω1=Inf)
    all_Ω = Float64[]
    all_ω1 = Float64[]
    all_t = Float64[]
    all_spectra = NMRData[]
    bf = 0.0

    for filename in filenames
        expt = loadnmr(filename)
        expt /= NMRTools.scale(expt)

        # Get field strength (first experiment sets it)
        if bf == 0.0
            bf = expt[1, :bf]  # Larmor frequency in MHz
        end

        # Extract conditions from annotations
        powers = annotations(expt, :r1rho, :power)
        offsets = annotations(expt, :r1rho, :offset)
        durations = annotations(expt, :r1rho, :duration)
        channel = annotations(expt, :r1rho, :channel)

        # get reference pulse
        p, pl = referencepulse(expt, channel)

        if isempty(all_spectra)
            error("No data loaded from provided files")
        end

        # Normalize spectra by maximum intensity
        mx = maximum(maximum(data(s)) for s in all_spectra)
        all_spectra = [s / mx for s in all_spectra]

        return LoadedData(all_Ω, all_ω1, all_t, all_spectra, bf)
    end
end

"""
    load_cest_data(filenames::Vector{String})

Load CEST experiments from NMR data files.

# Arguments
- `filenames`: Paths to Bruker experiment directories

# Returns
- `LoadedData` containing conditions and spectra
"""
function load_cest_data(filenames::Vector{String})
    all_Ω = Float64[]
    all_ω1 = Float64[]
    all_t = Float64[]
    all_spectra = NMRData[]
    bf = 0.0

    for filename in filenames
        expt = loadnmr(filename)
        expt /= NMRTools.scale(expt)

        # Get field strength
        if bf == 0.0
            bf = expt[1, :bf]
        end

        # Extract CEST conditions from annotations
        offsets = annotations(expt, :cest, :offset)
        power = annotations(expt, :cest, :power)
        sattime = annotations(expt, :cest, :saturation_time)
        channel = annotations(expt, :cest, :channel)

        # Convert power to rad/s
        p, pl = referencepulse(expt, channel)
        ν1_hz = hz(power, pl, p, 90)
        ω1 = 2π * ν1_hz

        # Convert offsets to rad/s
        Ω_hz = hz(offsets, dims(expt, F1Dim))
        Ω = 2π * Ω_hz

        # Build conditions
        nΩ = length(Ω)
        for i in 1:nΩ
            push!(all_Ω, Ω[i])
            push!(all_ω1, ω1)
            push!(all_t, sattime)
            push!(all_spectra, expt[:, i])
        end
    end

    if isempty(all_spectra)
        error("No CEST data loaded from provided files")
    end

    # Normalize
    mx = maximum(maximum(data(s)) for s in all_spectra)
    all_spectra = [s / mx for s in all_spectra]

    return LoadedData(all_Ω, all_ω1, all_t, all_spectra, bf)
end

"""
    integrate_spectra(loaded::LoadedData, roi, noiseroi)

Integrate loaded spectra over the specified region.

# Arguments
- `loaded`: LoadedData from load_r1rho_data or load_cest_data
- `roi`: Integration region (ppm interval, e.g., -120.5 .. -119.5)
- `noiseroi`: Noise region for uncertainty estimation

# Returns
- `intensities`: Integrated peak intensities
- `uncertainties`: Noise-based uncertainties
"""
function integrate_spectra(loaded::LoadedData, roi, noiseroi)
    intensities = Float64[]
    uncertainties = Float64[]

    for spec in loaded.spectra
        # Integrate signal region
        signal = sum(data(spec[roi]))
        push!(intensities, signal)

        # Estimate noise from noise region
        noise_data = data(spec[noiseroi])
        noise = std(noise_data) * sqrt(length(data(spec[roi])))
        push!(uncertainties, noise)
    end

    # Normalize
    mx = maximum(abs.(intensities))
    intensities ./= mx
    uncertainties ./= mx

    return intensities, uncertainties
end

"""
    make_r1rho_experiment(loaded::LoadedData, intensities, uncertainties)

Create an R1rhoExperiment from loaded data and integrated intensities.
"""
function make_r1rho_experiment(loaded::LoadedData, intensities::Vector{Float64},
                               uncertainties::Vector{Float64})
    return R1rhoExperiment(loaded.Ω, loaded.ω1, loaded.t, intensities, uncertainties,
                           loaded.bf)
end

"""
    make_cest_experiment(loaded::LoadedData, intensities, uncertainties)

Create a CESTExperiment from loaded data and integrated intensities.
"""
function make_cest_experiment(loaded::LoadedData, intensities::Vector{Float64},
                              uncertainties::Vector{Float64})
    return CESTExperiment(loaded.Ω, loaded.ω1, loaded.t, intensities, uncertainties,
                          loaded.bf)
end

"""
    filter_r1rho_files(filenames::Vector{String})

Filter filenames to only R1ρ experiments.
"""
function filter_r1rho_files(filenames::Vector{String})
    filter(filenames) do f
        try
            expt = loadnmr(f)
            hasannotations(expt) && "r1rho" in annotations(expt, :experiment_type)
        catch
            false
        end
    end
end

"""
    filter_cest_files(filenames::Vector{String})

Filter filenames to only CEST experiments.
"""
function filter_cest_files(filenames::Vector{String})
    filter(filenames) do f
        try
            expt = loadnmr(f)
            hasannotations(expt) && "cest" in annotations(expt, :experiment_type)
        catch
            false
        end
    end
end

"""
    filter_r1_calibration_files(filenames::Vector{String})

Filter filenames to R1 calibration experiments.
"""
function filter_r1_calibration_files(filenames::Vector{String})
    filter(filenames) do f
        try
            expt = loadnmr(f)
            hasannotations(expt) &&
                "relaxation" in annotations(expt, :experiment_type) &&
                "R1" in annotations(expt, :features)
        catch
            false
        end
    end
end
