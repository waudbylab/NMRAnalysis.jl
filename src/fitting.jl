# Fitting utilities shared by every module.

"""
    stderrors(fit) -> Vector{Float64}

Standard errors of a least-squares fit, or `NaN`s where the covariance matrix cannot be
inverted.

`LsqFit.stderror` inverts JᵀJ, which is singular whenever the data leave a parameter
unconstrained: an integration region dragged onto the baseline, an amplitude that has gone
to zero, two coincident points on the fit axis. Throwing is the right answer in a script
and the wrong one in an interactive fit, where the same exception arrives inside a Makie
callback and takes the window's refit down with it. A `NaN` uncertainty instead propagates
into the results and the plots, where it is visible and harmless.

The warning is capped: a fit re-run on every drag of a region would otherwise fill the
console with the same line.
"""
function stderrors(fit)
    try
        return stderror(fit)
    catch e
        e isa SingularException || e isa PosDefException || e isa LAPACKException ||
            rethrow()
        @warn "A fit's covariance matrix could not be inverted, so its uncertainties are " *
              "reported as NaN: the data do not constrain every parameter." maxlog=3
        return fill(NaN, length(coef(fit)))
    end
end
