"""
    viscosity(solvent, T)

Calculate the viscosity of a solvent (in mPa s) at a given temperature.

# Arguments
- `solvent`: Should be either `:h2o` or `:d2o`.
- `T`: The temperature (in K) at which the viscosity is calculated.

# Notes
The calculation is based on the formula proposed by Cho et al in the paper "J Phys Chem B (1999) 103 1991-1994".
"""
function viscosity(solvent, T)
    @info "Viscosity: calculation based on Cho et al, J Phys Chem B (1999) 103 1991-1994"
    if solvent==:h2o
        A = 802.25336
        a = 3.4741e-3
        b = -1.7413e-5
        c = 2.7719e-8
        γ = 1.53026
        T0 = 225.334
    elseif solvent==:d2o
        A = 885.60402
        a = 2.799e-3
        b = -1.6342e-5
        c = 2.9067e-8
        γ = 1.55255
        T0 = 231.832
    else
        @error "solvent not recognised (should be :h2o or :d2o)"
    end

    DT = T - T0

    return A * (DT + a*DT^2 + b*DT^3 + c*DT^4)^(-γ)
end
