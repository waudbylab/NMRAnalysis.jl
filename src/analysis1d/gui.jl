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

"""
    pickregion(traces; peakppm=nothing, noiseppm=nothing, ppmwidth=0.1) -> (; peakppm, noiseppm, ppmwidth)
    pickregion(specs; kwargs...)

Interactively choose a single integration region and a noise region, returning the
integration triple. All traces are overlaid, so a region can be chosen that suits every
spectrum at once - which is what is wanted when one region must serve a whole set of
experiments (e.g. `Exchange1D`'s `integrate!`, which applies one triple to every
experiment in the problem).

This is the region-selection front-end on its own, without any fitting: use it where the
heavy fitting lives elsewhere.
"""
function pickregion(traces::AbstractVector{Trace}; peakppm=nothing, noiseppm=nothing,
                    ppmwidth=0.1)
    GLMakie.activate!(; focus_on_show=true,
                      title="NMRAnalysis.jl: select integration region")
    cols = Makie.wong_colors()

    t1 = first(traces)
    peakppm = isnothing(peakppm) ? t1.δ[argmax(t1.y)] : peakppm
    noiseppm = isnothing(noiseppm) ? t1.δ[max(1, round(Int, 0.9 * length(t1.δ)))] : noiseppm

    pk = Observable(Float64(peakppm))
    nz = Observable(Float64(noiseppm))
    w = Observable(Float64(ppmwidth))

    fig = Figure(; size=(1000, 620))
    ax = Axis(fig[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title="Drag the shaded regions, then press Accept")
    hlines!(ax, [0]; color=:grey)
    lines!(ax, _overlay_points([Point2f.(t.δ, t.y) for t in traces]); color=(:grey, 0.5),
           label="All spectra")
    peakspan = vspan!(ax, lift((p, ww) -> p - ww / 2, pk, w),
                      lift((p, ww) -> p + ww / 2, pk, w); alpha=0.3, color=cols[2],
                      label="Peak")
    noisespan = vspan!(ax, lift((p, ww) -> p - ww / 2, nz, w),
                       lift((p, ww) -> p + ww / 2, nz, w); alpha=0.3, color=cols[4],
                       label="Noise")
    axislegend(ax; position=:lt)

    ctrl = fig[2, 1] = GridLayout()
    ctrl[1, 1] = Label(fig, "Integration width (ppm):")
    tb = ctrl[1, 2] = Textbox(fig; stored_string=string(round(w[]; digits=3)),
                              validator=Float64, width=100)
    on(s -> w[] = parse(Float64, s), tb.stored_string)
    accept = ctrl[1, 3] = Button(fig; label="Accept")
    done = Ref(false)
    on(_ -> done[] = true, accept.clicks)

    dragging = Ref(:nothing)
    on(events(ax).mousebutton) do ev
        ev.button == Mouse.left || return Consume(false)
        if ev.action == Mouse.press
            dragging[] = if mouseover(fig, peakspan)
                :peak
            elseif mouseover(fig, noisespan)
                :noise
            else
                :nothing
            end
            return Consume(dragging[] != :nothing)
        elseif ev.action == Mouse.release
            dragging[] = :nothing
            return Consume(false)
        end
    end
    on(events(fig).mouseposition; priority=2) do _
        if dragging[] == :peak
            pk[] = mouseposition(ax)[1]
            return Consume(true)
        elseif dragging[] == :noise
            nz[] = mouseposition(ax)[1]
            return Consume(true)
        end
        return Consume(false)
    end

    display(fig)
    while isopen(fig.scene) && !done[]
        sleep(0.1)
    end
    GLMakie.closeall()
    return (; peakppm=pk[], noiseppm=nz[], ppmwidth=w[])
end

function pickregion(specs::AbstractVector; kwargs...)
    return pickregion(reduce(vcat, traces_from_spec.(specs)); kwargs...)
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
