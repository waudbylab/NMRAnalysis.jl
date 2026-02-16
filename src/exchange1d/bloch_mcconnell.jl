# Bloch-McConnell equations for chemical exchange

"""
    build_bloch_mcconnell(model::AbstractModel, modelpars, spinpars, Ω, ω1)

Build the extended Bloch-McConnell Liouvillian matrix.

The magnetization vector is extended to handle the inhomogeneous R1 relaxation term:
    M = [Mx_1, My_1, Mz_1, ..., Mx_n, My_n, Mz_n, 1]^T

This converts the inhomogeneous ODE dM/dt = L*M + b into a homogeneous system
by augmenting with a constant element.

# Arguments
- `model`: Exchange model (NoExchange, TwoState, etc.)
- `modelpars`: Model parameters (kex, pB for TwoState)
- `spinpars`: Spin parameters (R1, R2, δ for each state)
- `Ω`: Carrier offset from spin frequency (rad/s)
- `ω1`: Spin-lock or saturation field strength (rad/s)

# Returns
- `L`: (3n+1) × (3n+1) Liouvillian matrix
"""
function build_bloch_mcconnell(model::AbstractModel, modelpars, spinpars, Ω, ω1)
    n = nstates(model)
    K = calculate_K(model, modelpars)
    p0 = calculate_p0(model, modelpars)
    R1 = value(spinpars.R1)

    # Extended Liouvillian: (3n+1) × (3n+1)
    L = zeros(3n + 1, 3n + 1)

    # Single-spin Bloch matrices for each state
    for i in 1:n
        idx = 3(i - 1) .+ (1:3)
        R2 = value(spinpars.R2[i])
        δ = value(spinpars.δ[i])

        # Effective offset in rotating frame
        Ω_eff = Ω - δ

        # Bloch matrix in rotating frame with spin-lock along x
        # dMx/dt = -R2*Mx - Ω_eff*My
        # dMy/dt = Ω_eff*Mx - R2*My - ω1*Mz
        # dMz/dt = ω1*My - R1*Mz + R1*Mz_eq
        L[idx, idx] = [-R2 -Ω_eff 0
                       Ω_eff -R2 -ω1
                       0 ω1 -R1]

        # Inhomogeneous term: R1 relaxation back to equilibrium
        # The equilibrium z-magnetization for state i is p0[i]
        L[3i, 3n + 1] = R1 * p0[i]
    end

    # The constant element (index 3n+1) doesn't evolve
    L[3n + 1, 3n + 1] = 0.0

    # Add exchange between corresponding magnetization components
    # K[i,j] is rate from j to i (gain term for i from j)
    # K[i,i] is negative of total rate out of i (loss term)
    for i in 1:n, j in 1:n
        for c in 1:3  # Mx, My, Mz components
            L[3(i - 1) + c, 3(j - 1) + c] += K[i, j]
        end
    end

    return L
end

"""
    simulate(model::AbstractModel, modelpars, spinpars, Ω, ω1, t, M0)

Simulate magnetization evolution under the Bloch-McConnell equations.

# Arguments
- `model`: Exchange model
- `modelpars`: Model parameters
- `spinpars`: Spin parameters
- `Ω`: Carrier offset (rad/s)
- `ω1`: Spin-lock field strength (rad/s)
- `t`: Evolution time (s)
- `M0`: Initial magnetization vector (length 3n)

# Returns
- `M`: Final magnetization vector (length 3n)
"""
function simulate(model::AbstractModel, modelpars, spinpars, Ω, ω1, t, M0)
    L = build_bloch_mcconnell(model, modelpars, spinpars, Ω, ω1)

    # Extend initial magnetization with constant element
    M0_ext = vcat(M0, 1.0)

    # Evolve using matrix exponential
    M_ext = exp(L * t) * M0_ext

    # Return magnetization components (drop the constant)
    return M_ext[1:(end - 1)]
end

"""
    observe_z(model::AbstractModel, M)

Extract the observable z-magnetization (sum over all states).

# Arguments
- `model`: Exchange model
- `M`: Magnetization vector (length 3n)

# Returns
- Total z-magnetization (scalar)
"""
function observe_z(model::AbstractModel, M)
    n = nstates(model)
    return sum(M[3i] for i in 1:n)
end

"""
    observe_effective(model::AbstractModel, M, Ω, ω1, spinpars; flip=false)

Extract the magnetization component along the effective field direction.

For R1ρ experiments, the detected signal is proportional to the magnetization
along the effective field, which is what actually decays exponentially.

# Arguments
- `model`: Exchange model
- `M`: Magnetization vector (length 3n)
- `Ω`: Carrier offset (rad/s)
- `ω1`: Spin-lock field strength (rad/s)
- `spinpars`: Spin parameters (for chemical shifts)
- `flip`: If true, observe along the flipped effective field (θ + 180°)

# Returns
- Total magnetization along effective field (scalar)
"""
function observe_effective(model::AbstractModel, M, Ω, ω1, spinpars; flip=false)
    n = nstates(model)
    total = 0.0

    for i in 1:n
        δ = value(spinpars.δ[i])
        Ω_eff = Ω - δ

        # Tilt angle of effective field from z-axis
        θ = atan(ω1, Ω_eff)

        # For flip=true, add 180° to observe along opposite direction
        if flip
            θ += π
        end

        # Magnetization components for this state
        Mx = M[3(i - 1) + 1]
        My = M[3(i - 1) + 2]
        Mz = M[3(i - 1) + 3]

        # Projection onto effective field direction
        # Effective field is in xz plane: (sin(θ), 0, cos(θ))
        M_eff = Mx * sin(θ) + Mz * cos(θ)
        total += M_eff
    end

    return total
end
