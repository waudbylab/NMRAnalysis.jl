# define a 'safe' square root function that returns 0 for negative inputs
safesqrt(x) = x < 0 ? NaN : sqrt(x)

function estimatekexfromK(K, σΔδ, bf)
    dw = 0 ± (σΔδ * 1e-6 * 2π * bf)
    kex = safesqrt(K^2 - dw^2)
    return kex = typeof(kex)(filter(!isnan, kex.particles)) # remove NaNs
end
