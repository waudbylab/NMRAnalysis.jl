"""
    plot_result(expt::R1Experiment, params; kwargs...)

Plot an R1 experiment with observed data (error bars), fitted curve, and residuals.

Upper panel: normalised intensities vs delay time with fit overlay.
Lower panel: weighted residuals (observed - predicted) / σ.
"""
function plot_result(expt::R1Experiment; kwargs...)
    tau = expt.delays
    obs = Measurements.value.(expt.observed_intensities)
    errs = Measurements.uncertainty.(expt.observed_intensities)
    pred = expt.predicted_intensities

    # fitted rate for title
    fl = field_label(expt)
    R1 = params.spin[Symbol("R1_", fl)][1]

    # smooth fitted curve
    x = LinRange(0, maximum(tau) * 1.1, 100)
    tag = Symbol("R1_", fl)
    I0 = params.nuisance[Symbol(tag, "_I0")]
    if expt.fitting_model == :inversion_recovery
        inv_factor = params.nuisance[Symbol(tag, "_inv_factor")]
        yfit = I0 .* (1 .- inv_factor .* exp.(-x .* R1))
    else
        yfit = I0 .* exp.(-x .* R1)
    end

    # weighted residuals
    wres = (obs .- pred) ./ errs

    # upper panel: data + fit
    p1 = scatter(tau, obs;
                 yerror=errs,
                 label="observed",
                 frame=:box,
                 xlabel="",
                 ylabel="Normalised intensity",
                 title="R1 = $(round(R1; digits=3)) s⁻¹",
                 grid=nothing,
                 kwargs...)
    plot!(p1, x, yfit;
          label="fit",
          z_order=:back)
    hline!(p1, [0]; color=:black, lw=0.5, primary=false)

    # lower panel: residuals
    p2 = scatter(tau, wres;
                 label=false,
                 frame=:box,
                 xlabel="Delay / s",
                 ylabel="Residual / σ",
                 grid=nothing,
                 markersize=4)
    hline!(p2, [0]; color=:black, lw=0.5, ls=:dash, primary=false)

    plt = plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]), link=:x)
    return plt
end

"""
    plot_result(prob::ExchangeProblem, fit_result)

Plot all experiments in the problem using the fitted parameters from `fit_result`
(as returned by `fit`). Returns a vector of plots, one per experiment.
"""
function plot_result(prob::ExchangeProblem, fit_result, kwargs...)
    params = fit_result.params_value
    simulate!(prob, params)  # update predicted_intensities for all experiments
    plots = [plot_result(expt; kwargs...) for expt in prob.experiments]
    return plots
end
