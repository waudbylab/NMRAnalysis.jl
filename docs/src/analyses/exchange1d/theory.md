# Theory and Calculation Methods

This page describes the Bloch-McConnell formalism and simulation methods shared by all
exchange models. For what each model represents physically and the population/binding
equations specific to it, see [Exchange Models](models.md).

## Bloch-McConnell Equations

Exchange analysis uses the Bloch-McConnell formalism, which extends the Bloch equations to multiple exchanging states. For ``N`` states, the magnetisation vector ``\mathbf{M}`` evolves as:

```math
\frac{d\mathbf{M}}{dt} = \mathbf{L} \, \mathbf{M}
```

where ``\mathbf{L}`` is the ``3N \times 3N`` Liouvillian superoperator incorporating:
- Chemical shift evolution (``\Omega_i`` for each state)
- Relaxation (``R_1``, ``R_{2,i}`` for each state)
- RF fields (spin-lock or saturation)
- Chemical exchange (kinetic rate matrix ``\mathbf{K}``)

## CEST Simulation

For CEST experiments, the Liouvillian is augmented with an inhomogeneous term to account for relaxation back to equilibrium. The predicted CEST profile is computed by propagating the initial equilibrium magnetisation through the saturation period:

```math
\mathbf{M}(T_\text{sat}) = \exp(\mathbf{L}_\text{inhom} \cdot T_\text{sat}) \, \mathbf{M}_0
```

The observed intensity is the sum of ``M_z`` components across all states.

## R1ρ Simulation

For on- and off-resonance R1ρ experiments, the effective relaxation rate at each spin-lock
condition is obtained directly from the Liouvillian ``\mathbf{L}`` (evaluated at that
experiment's spin-lock offset and field strength), via the inverse-trace relation

```math
R_{1\rho} = -\frac{1}{\mathrm{tr}(\mathbf{L}^{-1})},
```

sometimes referred to as the Koss method. On-resonance experiments vary the spin-lock field
strength at zero offset; off-resonance experiments hold the spin-lock field fixed and vary
the offset.

## R1 Simulation

R1 experiments are fitted independently of the exchange model (R1 decay is not sensitive to chemical exchange under typical conditions). Depending on the experiment type:

- **Exponential decay**: ``I(t) = I_0 \exp(-R_1 t)``
- **Inversion recovery**: ``I(t) = I_0 (1 - f \exp(-R_1 t))`` where ``f`` is the inversion factor

## Exchange Matrix

The kinetic exchange matrix ``\mathbf{K}`` encodes the rates of interconversion between states. Each element ``K_{ij}`` is the rate from state ``j`` to state ``i``, and column sums are zero (conservation of total magnetisation). ``\mathbf{K}`` is built once per model from the current parameters and (for binding models) the sample concentrations at each titration point — see [Exchange Models](models.md) for the matrix specific to each model.

## Binding Equilibria

The binding models additionally require the equilibrium populations and, where relevant, the free titrant concentration to be solved analytically from total concentrations at each titration point, rather than fitted directly. These are closed-form solutions of the coupled mass-action and mass-balance equations for each mechanism:

- **Two-state binding**: a single quadratic equation in the free titrant concentration.
- **Three-state binding**: a quadratic equation for two independent, non-competing 1:1 sites on the same molecule.
- **Induced fit binding**: a system of two equations (bimolecular association plus a first-order conformational step), reducing to a single quadratic in the same manner, ported from the TITAN `bmInducedFit` model.

See [Exchange Models](models.md) for the equations themselves. In each case, the resulting
species concentrations are normalised to their own sum (rather than divided through by the
nominal total concentration directly) so that population fractions remain exactly
normalised under floating-point rounding.
