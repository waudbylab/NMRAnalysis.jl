"""
    liouvillian(model, params, field_teslas, bf_hz, offset_hz, spinlock_hz, sampleconcentrations)

Build the 3N × 3N Bloch-McConnell Liouvillian matrix for N exchanging states.

The magnetisation vector is ordered as [Mx₁, My₁, Mz₁, ..., MxN, MyN, MzN].

# Arguments
- `model`: exchange model (determines number of states and kinetic matrix)
- `params`: ComponentArray with `model` and `spin` sections
- `field_teslas`: magnetic field strength in Tesla (for parameter lookup keys)
- `bf_hz`: base (Larmor) frequency in Hz for the observed nucleus (for ppm → Hz)
- `offset_hz`: spin-lock carrier offset in Hz
- `spinlock_hz`: spin-lock field amplitude in Hz
- `sampleconcentrations`: Dict mapping molecule names to concentrations (for binding models)

Chemical shifts are stored as ppm in `params.spin.delta` and converted to Hz
using `bf_hz`. R2 is per-state per-field in `params.spin.R2_<field>`.
R1 is shared across states per-field in `params.spin.R1_<field>` (length-1 vector).
"""
function liouvillian(model::AbstractModel, params, field_teslas::Float64, bf_hz::Float64,
                     offset_hz, spinlock_hz, sampleconcentrations)
    N = nstates(model)
    fl = field_label(field_teslas)

    L = zeros(3N, 3N)

    omega1 = 2π * spinlock_hz

    R2_key = Symbol("R2_", fl)
    R1_key = Symbol("R1_", fl)

    for i in 1:N
        # chemical shift offset in Hz, converted from ppm
        delta_ppm = params.spin.delta[i]
        delta_hz = 1e-6 * delta_ppm * bf_hz
        Omega_i = 2π * (delta_hz - offset_hz)

        # block indices for state i
        ix = 3(i - 1) + 1  # Mx
        iy = 3(i - 1) + 2  # My
        iz = 3(i - 1) + 3  # Mz

        # chemical shift evolution: couples Mx ↔ My
        L[ix, iy] = -Omega_i
        L[iy, ix] = Omega_i

        # spin-lock field along x: couples My ↔ Mz
        L[iy, iz] = -omega1
        L[iz, iy] = omega1

        # relaxation
        R2_i = params.spin[R2_key][i]
        R1_i = params.spin[R1_key][1]  # R1 is shared across states (stored as scalar)

        L[ix, ix] = -R2_i
        L[iy, iy] = -R2_i
        L[iz, iz] = -R1_i
    end

    # exchange: K ⊗ I₃
    K = exchange_matrix(model, params, sampleconcentrations)
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
    liouvsillian(model, params, expt::AbstractExperiment, offset_hz, spinlock_hz)

Convenience method that extracts field_teslas, bf_hz, and sample concentrations
from the experiment.
"""
function liouvillian(model::AbstractModel, params, expt::AbstractExperiment,
                     offset_hz, spinlock_hz)
    bf_hz = metadata(expt.spec, 1, :bf)
    return liouvillian(model, params, expt.field_teslas, bf_hz, offset_hz, spinlock_hz,
                       expt.sampleconcentrations)
end
