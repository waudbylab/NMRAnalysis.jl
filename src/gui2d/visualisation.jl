# Default trait - all experiments use cross sections unless specified
visualisationtype(::Experiment) = CrossSectionVisualisation()

# Plotting functions dispatch on the visualization trait
function plot_peak!(panel, peak, expt::Experiment)
    plot_peak!(panel, peak, expt, visualisationtype(expt))
end

function makepeakplot!(gui, state, expt::Experiment)
    makepeakplot!(gui, state, expt, visualisationtype(expt))
end

function completestate!(state, expt::Experiment)
    completestate!(state, expt, visualisationtype(expt))
end

function get_model_data(peak, expt::Experiment)
    error("get_model_data not implemented for $(typeof(expt))")
end

function save_peak_plots!(expt::E, folder::AbstractString) where E <: Experiment
    CairoMakie.activate!()
    
    for peak in expt.peaks[]
        fig = Figure()
        plot_peak!(fig, peak, expt)  # Uses the trait-based plot_peak! we defined
        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end
    
    GLMakie.activate!()
end


## Cross-section visualisation

# Per-plane multiplicative factor applied to displayed cross-section intensities. The default
# leaves them untouched; experiments whose stored `specdata.z` is not amplitude-comparable
# across planes override this (see MovingExperiment, which stores per-plane S/N).
crosssection_scale(::Experiment, i) = 1.0

function get_cross_section_data(peak, expt::Experiment)
    if isnothing(peak)
        return (Point2f[], Point2f[], Point2f[], Point2f[])
    end

    xobs, yobs, xfit, yfit = Vector{Point2f}[], Vector{Point2f}[], Vector{Point2f}[], Vector{Point2f}[]

    for i in 1:nslices(expt)
        x = expt.specdata.x[i]
        y = expt.specdata.y[i]
        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]

        # Get cross sections
        ix = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        iy = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        ix0 = findnearest(x, x0)
        iy0 = findnearest(y, y0)

        # Rescale obs and fit together so peak heights are comparable across planes (e.g.
        # differing numbers of scans / receiver gain); see crosssection_scale.
        s = crosssection_scale(expt, i)
        push!(xobs, Point2f.(x[ix], s .* expt.specdata.z[i][ix, iy0]))
        push!(yobs, Point2f.(y[iy], s .* expt.specdata.z[i][ix0, iy]))
        push!(xfit, Point2f.(x[ix], s .* expt.specdata.zfit[][i][ix, iy0]))
        push!(yfit, Point2f.(y[iy], s .* expt.specdata.zfit[][i][ix0, iy]))
    end
    
    return (flatten_with_nan_separator(xobs),
            flatten_with_nan_separator(yobs),
            flatten_with_nan_separator(xfit),
            flatten_with_nan_separator(yfit))
end

function makepeakplot!(gui, state, expt, ::CrossSectionVisualisation)
    @debug "making peak plot for cross section visualisation"
    gui[:axpeakplotX] = axX = Axis(gui[:panelpeakplot][1, 1];
                                   xlabel="δX / ppm",
                                   xreversed=true)
    gui[:axpeakplotY] = axY = Axis(gui[:panelpeakplot][1, 2];
                                   xlabel="δY / ppm",
                                   xreversed=true)

    hlines!(axX, [0]; linewidth=0)
    scatterlines!(axX, state[:peak_plot_xobs])
    lines!(axX, state[:peak_plot_xfit]; color=:red)
    hlines!(axY, [0]; linewidth=0)
    scatterlines!(axY, state[:peak_plot_yobs])
    lines!(axY, state[:peak_plot_yfit]; color=:red)
end

function plot_peak!(panel, peak, expt, ::CrossSectionVisualisation)
    xobs, yobs, xfit, yfit = get_cross_section_data(peak, expt)
    
    axX = Axis(panel[1, 1], xlabel="δX / ppm", xreversed=true)
    axY = Axis(panel[1, 2], xlabel="δY / ppm", xreversed=true)
    
    hlines!(axX, [0]; linewidth=0)
    scatterlines!(axX, xobs)
    lines!(axX, xfit; color=:red)
    hlines!(axY, [0]; linewidth=0)
    scatterlines!(axY, yobs)
    lines!(axY, yfit; color=:red)
end

function completestate!(state, expt, ::CrossSectionVisualisation)\
    @debug "completing state for cross section visualisation"
    state[:peak_plot_data] = lift(peak -> get_cross_section_data(peak, expt), state[:current_peak])
    state[:peak_plot_xobs] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_yobs] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_xfit] = lift(d -> d[3], state[:peak_plot_data])
    state[:peak_plot_yfit] = lift(d -> d[4], state[:peak_plot_data])
end


## Model fit visualisation

function plot_peak!(panel, peak, expt, ::ModelFitVisualisation)
    obs_points, obs_errors, fit_points, skip_points, skip_errors = get_model_data(peak, expt)

    ax = Axis(panel[1, 1],
              xlabel=get_model_xlabel(expt),
              ylabel=get_model_ylabel(expt))

    hlines!(ax, [0]; linewidth=0)
    lines!(ax, fit_points; color=:red)
    errorbars!(ax, obs_errors; whiskerwidth=10)
    scatter!(ax, obs_points)
    if !isempty(skip_points)
        errorbars!(ax, skip_errors; whiskerwidth=10, color=:gray60)
        scatter!(ax, skip_points; color=:transparent, strokecolor=:gray60,
                 strokewidth=1.5, markersize=8)
    end
end

function completestate!(state, expt, ::ModelFitVisualisation)
    @debug "completing state for model fit visualisation"
    state[:peak_plot_data] = lift(peak -> get_model_data(peak, expt), state[:current_peak])
    state[:peak_plot_obs]         = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_err]         = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_fit]         = lift(d -> d[3], state[:peak_plot_data])
    state[:peak_plot_skip_obs]    = lift(d -> d[4], state[:peak_plot_data])
    state[:peak_plot_skip_err]    = lift(d -> d[5], state[:peak_plot_data])
end

function makepeakplot!(gui, state, expt, ::ModelFitVisualisation)
    @debug "making peak plot for model fit visualisation"
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel=get_model_xlabel(expt),
                                 ylabel=get_model_ylabel(expt))

    hlines!(ax, [0]; linewidth=0)
    lines!(ax, state[:peak_plot_fit]; color=:red)
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    scatter!(ax, state[:peak_plot_obs])
    # Skipped planes: open grey markers
    errorbars!(ax, state[:peak_plot_skip_err]; whiskerwidth=10, color=:gray60)
    scatter!(ax, state[:peak_plot_skip_obs]; color=:transparent, strokecolor=:gray60,
             strokewidth=1.5, markersize=8)
end

get_model_xlabel(::Experiment) = "x"
get_model_ylabel(::Experiment) = "y"