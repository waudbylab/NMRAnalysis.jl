function model_R1rho_onres(νSL, p)
    R20 = p[2]
    Rex = p[3]
    if Rex < 0
        Rex = 0.0
    end
    lnK = p[4]
    @. R20 + Rex * exp(2lnK) / (exp(2lnK) + (2π * νSL)^2)
end
model_I_onres(TSL, νSL, p) = p[1] * exp.(-TSL .* model_R1rho_onres(νSL, p))
model_I_onres(x, p) = model_I_onres(x[:, 1], x[:, 2], p)

function fit_onres(state, p0)
    # p0 = [I0, R2_0, Rex, lnK]
    dataset = state[:dataset]

    x = hcat(dataset.TSLs, dataset.νSLs)
    intensities = state[:intensities][]

    @debug "running fit_onres..." p0
    fit = LsqFit.curve_fit(model_I_onres, x, intensities, p0; lower=[-Inf, 0.0, 0.0, -10.0],
                           upper=[Inf, 1000.0, 1000.0, 40.0])
    @debug "fit_onres results" fit.param fit.converged
    if !fit.converged
        @warn "fit did not converge"
    elseif any(isnan, fit.param)
        @warn "fit returned NaN parameters"
    end
    # if fit.param[3] < 0
    #     @warn "fit returned negative Rex; refitting with Rex=0 constraint"
    #     p0_constrained = copy(fit.param)
    #     p0_constrained[3] = 0.0
    #     fit = LsqFit.curve_fit(model_I_onres, x, intensities, p0_constrained; lower=[-Inf, -Inf, 0.0, -10.], upper=[Inf, 1000., 1000., 40.])
    #     @debug "refit results" fit.param fit.converged
    # end
    return fit
end

# Null model for R1rho (no exchange)
model_R1rho_onres_null(νSL, p) = @. p[2]  # Only R2,0 term, no Rex contribution
model_I_onres_null(TSL, νSL, p) = p[1] * exp.(-TSL .* model_R1rho_onres_null(νSL, p))
model_I_onres_null(x, p) = model_I_onres_null(x[:, 1], x[:, 2], p)

# p0 = [I0, R2,0]
function fit_onres_null(state, p0)
    dataset = state[:dataset]

    x = hcat(dataset.TSLs, dataset.νSLs)
    intensities = state[:intensities][]

    @debug "running null model fit..." p0
    # Only need parameters p[1] (I0) and p[2] (R2,0)
    fit = LsqFit.curve_fit(model_I_onres_null, x, intensities, p0)
    @debug "null model fit results" fit.param fit.converged
    return fit
end

function fitexp(t, I, p)
    I0 = p[1]
    R0 = 1 / mean(t) # initial guess based on average measured times
    # R0 = p[2] + p[3] / 2  # R2,0 + Rex/2
    # if R0 < 0
    #     R0 = p[2]
    # end
    p0 = [R0]

    @debug "running fitexp..." I0 R0

    model(t, p) = @. I0 * exp(-t * p[1])

    R, Re = try
        fit = LsqFit.curve_fit(model, t, I, p0)
        R = coef(fit)[1]
        Re = try
            stderror(fit)[1]
        catch
            0.0
        end
        R, Re
    catch
        R0, 0.0
    end

    @debug "fitexp results" R ± Re

    return R ± Re
end

function optimisewidth!(state)
    # make a range of dx values
    # set them, fit, and store the result
    # find the integration width giving minimum error in parameters
    # set the integration width to that value

    # check the range of peak maxima positions (in first increments) and ensure range covers this
    xs = find_maxima(state, state[:peakppm][], state[:dx][])
    # find the maximum and minimum x values
    min_x = minimum(xs)
    max_x = maximum(xs)
    x = (min_x + max_x) / 2
    state[:peakppm][] = x
    state[:gui][:text_peakppm].displayed_string[] = string(round(x; digits=3))

    dxmin = (max_x - min_x) * 2 # minimum dx = 2x range
    dx0 = state[:dx][]
    dxs = logrange(0.02, 50, 21) .* dx0
    dxmin = max(0.01, dxmin)
    dxs = filter(x -> dxmin <= x <= 5, dxs)

    # scoring function
    score(state) =
        try
            det(estimate_covar(state[:fit][]))
        catch e
            Inf
        end
    scores = map(dxs) do dx
        state[:dx][] = dx
        @debug dx, score(state)
        return score(state)
    end

    best = argmin(scores)
    @debug "best dx: $(dxs[best])"
    return state[:dx][] = dxs[best]
end

function find_maxima(state, start_x, dx)
    xs = []
    idx = map(x -> x[1], state[:series])
    for spec in state[:dataset].spectra[idx]
        x = find_maximum(spec, start_x, dx)
        push!(xs, x)
    end
    return xs
end

function find_maximum(spectrum, start_x, dx)
    x = data(spectrum, F1Dim)
    y = data(spectrum)

    # Validate inputs
    length(x) == length(y) || throw(ArgumentError("x and y lists must be the same length"))
    !isempty(x) || throw(ArgumentError("Input lists cannot be empty"))

    # Define the search range
    lower_bound = start_x - dx
    upper_bound = start_x + dx

    # Find indices of points within the search range
    range_indices = findall(xi -> lower_bound <= xi <= upper_bound, x)

    # Check if there are any points in the range
    if isempty(range_indices)
        throw(ArgumentError("No data points found within the specified range"))
    end

    # Find the maximum y value within the range
    max_y_index = argmax(y[range_indices])
    max_index = range_indices[max_y_index]

    # Return the x value at the maximum
    return x[max_index]
end
