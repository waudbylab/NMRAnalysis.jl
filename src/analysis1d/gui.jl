# Interactive single-window GUI for the 1D experiments. Mirrors the idioms of the R1ρ and
# GUI2D GUIs (GLMakie figure, draggable shaded regions, plane slider, live `lift` refit).
# A second (exchange) window is deliberately out of scope here — none of relaxation /
# TRACT / nutation / STD / kinetics needs it.

"""
    gui!(expt::Experiment1D)

Launch the interactive analysis window for `expt`. The left column overlays all spectral
planes with draggable integration region(s) and a noise region; the result panel below
shows the live fit for the active region. Returns the GUI state when the window closes.
"""
function gui!(expt::Experiment1D)
    GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: 1D analysis")
    state = prepare_state(expt)
    state[:gui] = Dict{Symbol,Any}()
    gui = state[:gui]
    cols = Makie.wong_colors()

    fig = Figure(; size=(1200, 800))
    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    # --- spectral overlay ---
    ax = Axis(left[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title="Spectra – drag shaded regions to reposition")
    gui[:ax_spec] = ax
    hlines!(ax, [0]; color=:grey)
    lines!(ax, state[:overlay]; color=(:grey, 0.35), label="All planes")
    lines!(ax, state[:currentspectrum]; color=cols[1], label="Current plane")

    gui[:regionspans] = Any[]
    for i in 1:length(state[:regions][])
        lo = lift(rs -> rs[i].c - rs[i].w / 2, state[:regions])
        hi = lift(rs -> rs[i].c + rs[i].w / 2, state[:regions])
        push!(gui[:regionspans], vspan!(ax, lo, hi; alpha=0.3, color=cols[2]))
    end
    nlo = lift((c, w) -> c - w / 2, state[:noisec], state[:noisew])
    nhi = lift((c, w) -> c + w / 2, state[:noisec], state[:noisew])
    gui[:noisespan] = vspan!(ax, nlo, nhi; alpha=0.3, color=cols[4])
    axislegend(ax; position=:lt)

    # --- result panel ---
    xl, yl = result_labels(expt)
    ax2 = Axis(left[2, 1]; xlabel=xl, ylabel=yl, title="Fit (active region)")
    gui[:ax_fit] = ax2
    hlines!(ax2, [0]; linewidth=0)
    errorbars!(ax2, state[:ploterrors]; whiskerwidth=8)
    scatter!(ax2, state[:plotpoints]; color=cols[1], label="Observed")
    lines!(ax2, state[:plotfit]; color=cols[2], label="Fit")
    rowsize!(left, 1, Relative(0.55))

    # --- controls ---
    r = 0
    r += 1
    right[r, 1] = Label(fig, "Plane:")
    sl = right[r, 2] = Slider(fig; range=1:state[:nplanes], width=170)
    connect!(state[:currentplane], sl.value)

    r += 1
    right[r, 1] = Label(fig, "Region:")
    menu = right[r, 2] = Menu(fig; options=[r.label for r in state[:regions][]], width=170)
    on(menu.i_selected) do i
        return state[:active][] = i
    end

    r += 1
    right[r, 1] = Label(fig, "Region width (ppm):")
    tw = right[r, 2] = Textbox(fig; validator=Float64, width=170,
                               stored_string=string(round(state[:regions][][1].w; digits=3)))
    on(tw.stored_string) do s
        return set_region_width!(state, state[:active][], parse(Float64, s))
    end
    on(state[:active]) do i
        return tw.displayed_string[] = string(round(state[:regions][][i].w; digits=3))
    end

    r += 1
    right[r, 1] = Label(fig, "Noise width (ppm):")
    tnw = right[r, 2] = Textbox(fig; validator=Float64, width=170,
                                stored_string=string(round(state[:noisew][]; digits=3)))
    on(tnw.stored_string) do s
        return state[:noisew][] = parse(Float64, s)
    end

    r += 1
    btn_fit = right[r, 1:2] = Button(fig; label="Fitting: on")
    on(btn_fit.clicks) do _
        state[:isfitting][] = !state[:isfitting][]
        return btn_fit.label[] = state[:isfitting][] ? "Fitting: on" : "Fitting: off"
    end

    r += 1
    right[r, 1:2] = Label(fig, state[:summary]; tellwidth=false, halign=:left,
                          justification=:left)

    r += 1
    right[r, 1] = Label(fig, "Output folder:")
    tout = right[r, 2] = Textbox(fig; stored_string="out", width=170)
    on(tout.stored_string) do s
        return state[:outputdir][] = s
    end
    r += 1
    btn_save = right[r, 1:2] = Button(fig; label="Save results")
    on(btn_save.clicks) do _
        return save_results(state)
    end

    setup_dragging!(fig, ax, state)

    display(fig)
    autolimits!(ax2)
    while isopen(fig.scene)
        sleep(0.1)
    end
    GLMakie.closeall()
    return state
end

"""Wire mouse dragging of the region/noise spans on the spectrum axis."""
function setup_dragging!(fig, ax, state)
    gui = state[:gui]
    gui[:dragging] = :nothing
    on(events(ax).mousebutton) do ev
        ev.button == Mouse.left || return Consume(false)
        if ev.action == Mouse.press
            gui[:dragging] = :nothing
            for (i, sp) in enumerate(gui[:regionspans])
                if mouseover(fig, sp)
                    gui[:dragging] = (:region, i)
                    state[:active][] = i
                    break
                end
            end
            if gui[:dragging] == :nothing && mouseover(fig, gui[:noisespan])
                gui[:dragging] = :noise
            end
            return Consume(gui[:dragging] != :nothing)
        elseif ev.action == Mouse.release
            gui[:dragging] = :nothing
            return Consume(false)
        end
    end
    on(events(fig).mouseposition; priority=2) do _
        d = gui[:dragging]
        if d == :noise
            state[:noisec][] = mouseposition(ax)[1]
            return Consume(true)
        elseif d isa Tuple && d[1] == :region
            set_region_center!(state, d[2], mouseposition(ax)[1])
            return Consume(true)
        end
        return Consume(false)
    end
end

"""Save the active-region fit figure and a text summary to the output folder."""
function save_results(state)
    dir = joinpath(pwd(), state[:outputdir][])
    isdir(dir) || mkpath(dir)
    expt = state[:expt]
    cols = Makie.wong_colors()

    points, errors, fitline = state[:plotdata][]
    xl, yl = result_labels(expt)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=xl, ylabel=yl)
    hlines!(ax, [0]; linewidth=0)
    errorbars!(ax, errors; whiskerwidth=8)
    scatter!(ax, points; color=cols[1], label="Observed")
    lines!(ax, fitline; color=cols[2], label="Fit")
    axislegend(ax; position=:rt)
    save(joinpath(dir, "fit.pdf"), fig; backend=CairoMakie)

    open(joinpath(dir, "summary.txt"), "w") do f
        return println(f, state[:summary][])
    end
    @info "Saved results to $dir"
    return dir
end
