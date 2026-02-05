# Plotting functions for exchange analysis

using Plots

"""
    plot_r1rho_intensities(exp::R1rhoExperiment, result::FitResult)

Plot R1ρ intensity decay curves with fitted predictions.
"""
function plot_r1rho_intensities(exp::R1rhoExperiment, result::FitResult)
    # Group data by (Ω, ω1) combinations
    conditions = unique(zip(exp.Ω, exp.ω1))

    p = plot(xlabel="Relaxation time / s",
             ylabel="Intensity (normalized)",
             title="R1ρ Intensity Decay",
             legend=:topright,
             frame=:box)

    pred = result.predictions[1]  # Assuming single experiment

    for (Ω, ω1) in conditions
        idx = findall((exp.Ω .== Ω) .& (exp.ω1 .== ω1))
        t = exp.t[idx]
        I = exp.intensity[idx]
        σ = exp.σ[idx]
        I_pred = pred[idx]

        # Sort by time
        order = sortperm(t)
        t = t[order]
        I = I[order]
        σ = σ[order]
        I_pred = I_pred[order]

        # Label
        ω1_hz = round(ω1 / 2π, digits=0)
        Ω_hz = round(Ω / 2π, digits=0)
        label = Ω ≈ 0 ? "ω₁ = $ω1_hz Hz" : "Ω = $Ω_hz Hz, ω₁ = $ω1_hz Hz"

        # Plot observed with error bars
        scatter!(p, t, I, yerror=σ, label=label, markersize=4)

        # Plot predicted
        t_fine = LinRange(0, maximum(t) * 1.1, 50)
        # For smooth curve, we'd need to interpolate predictions
        plot!(p, t, I_pred, label="", linestyle=:solid)
    end

    return p
end

"""
    plot_cest_profile(exp::CESTExperiment, result::FitResult)

Plot CEST profile with fitted predictions.
"""
function plot_cest_profile(exp::CESTExperiment, result::FitResult)
    # Convert to ppm for display (assuming known B0 and nucleus)
    Ω_hz = exp.Ω ./ 2π

    pred = result.predictions[1]

    # Sort by offset
    order = sortperm(Ω_hz)

    p = plot(xlabel="Saturation offset / Hz",
             ylabel="Intensity (normalized)",
             title="CEST Profile",
             legend=:topright,
             frame=:box)

    scatter!(p, Ω_hz[order], exp.intensity[order],
             yerror=exp.σ[order],
             label="Observed",
             markersize=4)

    plot!(p, Ω_hz[order], pred[order],
          label="Fit",
          linewidth=2)

    return p
end

"""
    plot_residuals(exp::AbstractExperiment, result::FitResult)

Plot fit residuals.
"""
function plot_residuals(exp::AbstractExperiment, result::FitResult)
    pred = result.predictions[1]
    resid = (exp.intensity .- pred) ./ exp.σ

    p = plot(xlabel="Data point",
             ylabel="Residual (σ)",
             title="Fit Residuals",
             legend=false,
             frame=:box)

    scatter!(p, 1:length(resid), resid, markersize=4)
    hline!(p, [0], color=:black, linestyle=:dash)
    hline!(p, [-2, 2], color=:red, linestyle=:dot, alpha=0.5)

    return p
end

"""
    plot_results(experiments::Vector{<:AbstractExperiment}, results::Dict{<:AbstractModel, FitResult})

Create summary plots for all fitted models.
"""
function plot_results(experiments::Vector{<:AbstractExperiment},
                      results::Dict{<:AbstractModel, FitResult})
    plots = []

    # Plot for best model (by χ²)
    best_model = argmin(r -> r.χ2, values(results))
    best_result = results[best_model]

    for (i, exp) in enumerate(experiments)
        if exp isa R1rhoExperiment
            # Create a single-experiment result for plotting
            single_result = FitResult(
                best_result.model,
                best_result.modelpars,
                best_result.spinpars,
                best_result.χ2,
                best_result.ndata,
                best_result.npars,
                [best_result.predictions[i]]
            )
            push!(plots, plot_r1rho_intensities(exp, single_result))
        elseif exp isa CESTExperiment
            single_result = FitResult(
                best_result.model,
                best_result.modelpars,
                best_result.spinpars,
                best_result.χ2,
                best_result.ndata,
                best_result.npars,
                [best_result.predictions[i]]
            )
            push!(plots, plot_cest_profile(exp, single_result))
        end
    end

    # Add residual plot
    for (i, exp) in enumerate(experiments)
        single_result = FitResult(
            best_result.model,
            best_result.modelpars,
            best_result.spinpars,
            best_result.χ2,
            best_result.ndata,
            best_result.npars,
            [best_result.predictions[i]]
        )
        push!(plots, plot_residuals(exp, single_result))
    end

    # Combine plots
    n = length(plots)
    if n == 1
        return plots[1]
    elseif n == 2
        return plot(plots..., layout=(1, 2), size=(900, 400))
    else
        ncol = min(n, 2)
        nrow = ceil(Int, n / ncol)
        return plot(plots..., layout=(nrow, ncol), size=(450 * ncol, 350 * nrow))
    end
end
