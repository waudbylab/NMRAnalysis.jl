"""
    liouvillian(model, params, expt::AbstractExperiment, offset_hz, spinlock_hz)

Form the Liouvillian matrix for a given model, parameters, and experiment.
"""
function liouvillian(model, params, expt, spinlock_ppm, spinlock_hz)
    fl = field_label(expt)

    Δδ = params.spin.delta .- spinlock_ppm
    Ωi = Δδ .* expt.spec[1, :bf] * 1e-6
    ω = 2π * Ωi

    R1_key = Symbol("R1_", fl)
    R2_key = Symbol("R2_", fl)
    R1 = params.spin[R1_key][1]
    R2 = params.spin[R2_key]

    ω1 = 2π * spinlock_hz

    # @info Δδ, Ωi, ω, R1, R2, ω1

    N = nstates(model)
    L = zeros(3N, 3N)

    for i in 1:N
        # block indices for state i
        ix = 3(i - 1) + 1  # Mx
        iy = 3(i - 1) + 2  # My
        iz = 3(i - 1) + 3  # Mz

        # chemical shift evolution: couples Mx ↔ My
        L[ix, iy] = ω[i]
        L[iy, ix] = -ω[i]

        # spin-lock field along x: couples My ↔ Mz
        L[iy, iz] = -ω1
        L[iz, iy] = ω1

        # relaxation
        L[ix, ix] = -R2[i]
        L[iy, iy] = -R2[i]
        L[iz, iz] = -R1
    end

    # exchange: K ⊗ I₃
    K = exchangematrix(model, params, expt)
    for i in 1:N, j in 1:N
        kij = K[i, j]
        if kij != 0
            for d in 0:2  # Mx, My, Mz all exchange identically
                L[3(i - 1) + 1 + d, 3(j - 1) + 1 + d] += kij
            end
        end
    end

    return L
end

"""
    liouvillian_inhom(model, params, expt::AbstractExperiment, offset_hz, spinlock_hz)

Form the inhomogeneous Liouvillian matrix for a given model, parameters, and experiment.
"""
function liouvillian_inhom(model, params, expt, spinlock_ppm, spinlock_hz)
    fl = field_label(expt)

    Δδ = params.spin.delta .- spinlock_ppm
    Ωi = Δδ .* expt.spec[1, :bf] * 1e-6
    ω = 2π * Ωi

    R1_key = Symbol("R1_", fl)
    R2_key = Symbol("R2_", fl)
    R1 = params.spin[R1_key][1]
    R2 = params.spin[R2_key]

    ω1 = 2π * spinlock_hz

    # @info Δδ, Ωi, ω, R1, R2, ω1

    N = nstates(model)
    L = zeros(3N + 1, 3N + 1)

    for i in 1:N
        # block indices for state i
        ix = 3(i - 1) + 1  # Mx
        iy = 3(i - 1) + 2  # My
        iz = 3(i - 1) + 3  # Mz

        # chemical shift evolution: couples Mx ↔ My
        L[ix, iy] = ω[i]
        L[iy, ix] = -ω[i]

        # spin-lock field along x: couples My ↔ Mz
        L[iy, iz] = -ω1
        L[iz, iy] = ω1

        # relaxation
        L[ix, ix] = -R2[i]
        L[iy, iy] = -R2[i]
        L[iz, iz] = -R1
    end

    # exchange: K ⊗ I₃
    K = exchangematrix(model, params, expt)
    for i in 1:N, j in 1:N
        kij = K[j, i]
        if kij != 0
            for d in 0:2  # Mx, My, Mz all exchange identically
                L[3(i - 1) + 1 + d, 3(j - 1) + 1 + d] += kij
            end
        end
    end

    p0 = populations(model, params, expt)
    for i in 1:N
        L[3(i - 1) + 3, 3N + 1] = R1 * p0[i]  # Mz initialised to population
    end
    # @info L, K, p0
    return L
end
