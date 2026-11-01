# How chemical shifts are referenced: to a multiplet picked in the same analysis, or
# indirectly, through the IUPAC unified scale, to DSS picked in a ¹H analysis.

"""
    Referencing

How the shifts of an analysis are referenced. See [`PickedReference`](@ref) and
[`IndirectReference`](@ref).
"""
abstract type Referencing end

"""
    PickedReference(shift)

Shifts referenced to the centre of the reference multiplet picked in the same analysis,
which is placed at `shift` (ppm) in every plane.
"""
struct PickedReference <: Referencing
    shift::Float64
end

"""
    IndirectReference(frequencies, label)

Shifts on the IUPAC unified scale, referenced indirectly to DSS picked in a ¹H analysis of
the same samples. `frequencies` holds the frequency of 0 ppm on the ¹H scale in each plane
(Hz, `NaN` where DSS was not picked), and `label` names the ¹H reference.

On the unified scale the frequency of 0 ppm for nucleus X is Ξ_X times that of ¹H in DSS,
so a line at frequency ν has

    δ = (ν / (Ξ_X ν₀(¹H)) − 1) × 10⁶

where ν = bf (1 + δ′ × 10⁻⁶) for a line at apparent shift δ′ on an axis whose 0 ppm is at
`bf`. Because each spectrum's own `bf` enters, the result does not depend on how the X
spectrum was referenced when it was processed.

Harris, Becker, Cabral de Menezes, Goodfellow and Granger, Pure Appl. Chem. 73, 1795
(2001); Harris et al., Pure Appl. Chem. 80, 59 (2008).
"""
struct IndirectReference <: Referencing
    frequencies::Vector{Float64}
    label::String
end

"Largest relative difference between Ξ × ¹H base frequency and an X base frequency taken as the same field."
const FIELDTOLERANCE = 1e-4
