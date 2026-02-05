# Forward model predictions for exchange experiments

"""
    initial_magnetisation_r1rho(model::AbstractModel, modelpars, spinpars, Ω, ω1)

Calculate initial magnetization for R1ρ experiment.

The magnetization is aligned along the effective field in each state,
which is the standard preparation for R1ρ experiments.

# Arguments
- `model`: Exchange model
- `modelpars`: Model parameters
- `spinpars`: Spin parameters
- `Ω`: Carrier offset (rad/s)
- `ω1`: Spin-lock field strength (rad/s)

# Returns
- `M0`: Initial magnetization vector (length 3n)
"""
function initial_magnetisation_r1rho(model::AbstractModel, modelpars, spinpars, Ω, ω1)
    n = nstates(model)
    p0 = calculate_p0(model, modelpars)

    M0 = zeros(3n)
    for i in 1:n
        δ = value(spinpars.δ[i])
        Ω_eff = Ω - δ

        # Tilt angle of effective field from z-axis
        θ = atan(ω1, Ω_eff)

        # Magnetization aligned along effective field
        # In the tilted frame: M_eff = p0[i]
        # In the lab frame: Mx = p0[i]*sin(θ), Mz = p0[i]*cos(θ)
        M0[3(i - 1) + 1] = p0[i] * sin(θ)  # Mx
        M0[3(i - 1) + 2] = 0.0              # My
        M0[3(i - 1) + 3] = p0[i] * cos(θ)  # Mz
    end

    return M0
end

"""
    initial_magnetisation_cest(model::AbstractModel, modelpars, spinpars)

Calculate initial magnetization for CEST experiment.

The magnetization starts at thermal equilibrium (z-magnetization only).

# Arguments
- `model`: Exchange model
- `modelpars`: Model parameters
- `spinpars`: Spin parameters

# Returns
- `M0`: Initial magnetization vector (length 3n)
"""
function initial_magnetisation_cest(model::AbstractModel, modelpars, spinpars)
    n = nstates(model)
    p0 = calculate_p0(model, modelpars)

    M0 = zeros(3n)
    for i in 1:n
        M0[3i] = p0[i]  # Equilibrium z-magnetization
    end

    return M0
end

"""
    predict(exp::R1rhoExperiment, model::AbstractModel, modelpars, spinpars)

Compute predicted intensities for an R1ρ experiment.

Works whether the experiment has observed data or not (simulation mode).

The detected signal is proportional to the magnetization component along the
effective field direction, which decays exponentially with rate R1ρ.

# Arguments
- `exp`: R1ρ experiment with conditions
- `model`: Exchange model
- `modelpars`: Model parameters (kex, pB for TwoState)
- `spinpars`: Spin parameters (R1, R2, δ, amplitude)

# Returns
- Vector of predicted intensities
"""
function predict(exp::R1rhoExperiment, model::AbstractModel, modelpars, spinpars)
    amp = value(spinpars.R1rho_amplitude)

    map(exp.Ω, exp.ω1, exp.t) do Ω, ω1, t
        M0 = initial_magnetisation_r1rho(model, modelpars, spinpars, Ω, ω1)
        M = simulate(model, modelpars, spinpars, Ω, ω1, t, M0)
        amp * observe_effective(model, M, Ω, ω1, spinpars)
    end
end

"""
    predict(exp::CESTExperiment, model::AbstractModel, modelpars, spinpars)

Compute predicted intensities for a CEST experiment.

Works whether the experiment has observed data or not (simulation mode).

# Arguments
- `exp`: CEST experiment with conditions
- `model`: Exchange model
- `modelpars`: Model parameters (kex, pB for TwoState)
- `spinpars`: Spin parameters (R1, R2, δ, amplitude)

# Returns
- Vector of predicted intensities
"""
function predict(exp::CESTExperiment, model::AbstractModel, modelpars, spinpars)
    amp = value(spinpars.CEST_amplitude)

    map(exp.Ω, exp.ω1, exp.t) do Ω, ω1, t
        M0 = initial_magnetisation_cest(model, modelpars, spinpars)
        M = simulate(model, modelpars, spinpars, Ω, ω1, t, M0)
        amp * observe_z(model, M)
    end
end
