# Data loading from NMR experiments

using NMRTools
using Statistics: std

"""
    LoadedData

Container for loaded experiment data before integration.
"""
struct LoadedData
    Ω::Vector{Float64}      # Offsets (rad/s)
    ω1::Vector{Float64}     # Field strengths (rad/s)
    t::Vector{Float64}      # Times (s)
    spectra::Vector{NMRData}  # 1D spectra
    B0::Float64             # Field strength (T)
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
    B0 = 0.0

    for filename in filenames
        expt = loadnmr(filename)
        expt /= NMRTools.scale(expt)

        # Get field strength from first experiment
        if B0 == 0.0
            # B0 in Tesla from Larmor frequency
            # bf is in MHz, γ for 1H is 42.577 MHz/T
            bf = expt[1, :bf]  # MHz
            B0 = bf / 42.577  # Approximate, adjust for actual nucleus
        end

        # Extract conditions from annotations
        powers = annotations(expt, :r1rho, :power)
        offsets = annotations(expt, :r1rho, :offset)
        durations = annotations(expt, :r1rho, :duration)
        channel = annotations(expt, :r1rho, :channel)

        # Convert power to rad/s using reference pulse
        p, pl = referencepulse(expt, channel)
        ν1_hz = hz.(powers, pl, p, 90)  # Hz
        ω1 = 2π .* ν1_hz  # rad/s

        # Convert offsets to rad/s
        # offsets are typically in Hz or ppm, need to check format
        if offsets isa Number
            # Single offset (on-resonance)
            Ω_hz = [Float64(offsets)]
        else
            # Multiple offsets - convert from ppm or Hz
            Ω_hz = hz.(offsets, dims(expt, F1Dim))
        end
        Ω = 2π .* Ω_hz  # rad/s

        # Build grid of conditions
        nΩ = length(Ω)
        nω1 = length(ω1)
        nt = length(durations)

        # Determine array organization based on experiment dimensions
        if ndims(expt) == 3
            # 3D: (F1, ω1/Ω, t) or similar
            for i in 1:nt, j in 1:nΩ
                # Filter by spin-lock strength
                ω1_val = nω1 == 1 ? ω1[1] : ω1[j]
                if minω1 <= ω1_val <= maxω1
                    push!(all_Ω, Ω[j])
                    push!(all_ω1, ω1_val)
                    push!(all_t, durations[i])
                    push!(all_spectra, expt[:, j, i])
                end
            end
        elseif ndims(expt) == 2
            # 2D: (F1, combined index)
            for j in 1:nω1, i in 1:nt
                if minω1 <= ω1[j] <= maxω1
                    Ω_val = nΩ == 1 ? Ω[1] : 0.0  # On-resonance if single offset
                    push!(all_Ω, Ω_val)
                    push!(all_ω1, ω1[j])
                    push!(all_t, durations[i])
                    push!(all_spectra, expt[:, (j - 1) * nt + i])
                end
            end
        end
    end

    if isempty(all_spectra)
        error("No data loaded from provided files")
    end

    # Normalize spectra by maximum intensity
    mx = maximum(maximum(data(s)) for s in all_spectra)
    all_spectra = [s / mx for s in all_spectra]

    return LoadedData(all_Ω, all_ω1, all_t, all_spectra, B0)
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
    B0 = 0.0

    for filename in filenames
        expt = loadnmr(filename)
        expt /= NMRTools.scale(expt)

        # Get field strength
        if B0 == 0.0
            bf = expt[1, :bf]
            B0 = bf / 42.577
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
        Ω_hz = hz.(offsets, dims(expt, F1Dim))
        Ω = 2π .* Ω_hz

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

    return LoadedData(all_Ω, all_ω1, all_t, all_spectra, B0)
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
    R1rhoExperiment(loaded.Ω, loaded.ω1, loaded.t, intensities, uncertainties, loaded.B0)
end

"""
    make_cest_experiment(loaded::LoadedData, intensities, uncertainties)

Create a CESTExperiment from loaded data and integrated intensities.
"""
function make_cest_experiment(loaded::LoadedData, intensities::Vector{Float64},
                               uncertainties::Vector{Float64})
    CESTExperiment(loaded.Ω, loaded.ω1, loaded.t, intensities, uncertainties, loaded.B0)
end

"""
    filter_r1rho_files(filenames::Vector{String})

Filter filenames to only R1ρ experiments.
"""
function filter_r1rho_files(filenames::Vector{String})
    filter(filenames) do f
        try
            expt = loadnmr(f)
            hasannotations(expt) && "r1rho" in annotations(expt, :types)
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
            hasannotations(expt) && "cest" in annotations(expt, :types)
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
                "relaxation" in annotations(expt, :types) &&
                "R1" in annotations(expt, :features)
        catch
            false
        end
    end
end
