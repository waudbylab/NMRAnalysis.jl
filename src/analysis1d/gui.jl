# Interactive single-window GUI for the 1D experiments. Mirrors the idioms of the R1ρ and
# GUI2D GUIs (GLMakie figure, draggable shaded regions, spectrum slider, live `lift`
# refit, keyboard add/rename/delete matching the 2D peak-picking interface). A second
# (exchange) window is deliberately out of scope here — none of relaxation / TRACT /
# nutation / diffusion / STD / kinetics needs it.

# A tap vs. a drag when adding a region with 'A' - shorter moves than this (ppm) are
# treated as a single click, giving the default-width region.
const ADD_TAP_THRESHOLD = 0.005

"""
    gui!(expt::Experiment1D)

Launch the interactive analysis window for `expt`. The left column overlays all spectral
planes with draggable integration region(s) and a noise marker; the result panel below
shows the live fit for the active region. Returns the GUI state when the window closes.

# Keyboard shortcuts (mirroring the 2D fitting GUI)
- `A`: add a region. A quick tap adds a default-width (0.05 ppm) region under the
  cursor; press, drag, and release to choose the width.
- `R`: rename the active region.
- `D`: delete the active region.
- Left/Right arrows: step through spectra.
"""
function gui!(expt::Experiment1D)
    GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: 1D analysis")
    state = prepare_state(expt)
    state[:gui] = Dict{Symbol,Any}()
    gui = state[:gui]
    cols = Makie.wong_colors()
    singleregion = length(state[:regions][]) ≤ 1

    fig = Figure(; size=(1200, 800))
    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    # Without explicit column/row proportions, Makie sizes each cell to fit its content
    # and leaves the remainder of the window blank on resize - match GUI2D's layout and
    # give both columns a stretchy relative share so the interface fills the window.
    colsize!(fig.layout, 1, Auto(false, 3))
    colsize!(fig.layout, 2, Auto(false, 1))

    # mode-based background tint, matching the 2D fitting GUI's feedback style
    fig.scene.backgroundcolor = lift(state[:mode]) do mode
        if mode == :addingdrag
            RGBAf(1.0, 0.95, 0.7, 1.0)   # light yellow
        elseif mode == :renaming || mode == :renamingstart
            RGBAf(0.75, 0.94, 1.0, 1.0)  # light blue
        else
            RGBAf(1.0, 1.0, 1.0, 1.0)
        end
    end

    # --- spectral overlay ---
    ax = Axis(left[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title="Spectra – press A to add a region, drag to reposition")
    gui[:ax_spec] = ax
    hlines!(ax, [0]; color=:grey)
    lines!(ax, state[:overlay]; color=(:grey, 0.35), label="All spectra")
    lines!(ax, state[:currentspectrum]; color=cols[1], label="Current spectrum")

    gui[:regionspans] = Any[]
    rebuild_regionspans!(ax, state)
    on(rs -> length(rs) != length(gui[:regionspans]) && rebuild_regionspans!(ax, state),
       state[:regions])

    text!(ax, lift(rs -> [Point2f(reg.c, 0.0) for reg in rs], state[:regions]);
          text=lift(rs -> [reg.label for reg in rs], state[:regions]),
          align=(:center, :bottom), fontsize=11, offset=(0, 4))

    # noise marker: a vline for visibility, plus a narrow, mostly-transparent band as the
    # drag target (a true zero-width line is fiddly to grab with the mouse)
    vlines!(ax, state[:noisec]; color=:orchid, linewidth=2, label="Noise")
    gui[:noisehit] = vspan!(ax, lift(c -> c - 0.01, state[:noisec]),
                            lift(c -> c + 0.01, state[:noisec]); color=NOISE_COLOR)
    axislegend(ax; position=:lt)

    # add-region drag preview
    addlo = lift((m, a, c) -> m == :addingdrag ? min(a, c) : 0.0, state[:mode],
                 state[:addstart], state[:addcurrent])
    addhi = lift((m, a, c) -> m == :addingdrag ? max(a, c) : 0.0, state[:mode],
                 state[:addstart], state[:addcurrent])
    addalpha = lift(m -> m == :addingdrag ? 0.25 : 0.0, state[:mode])
    vspan!(ax, addlo, addhi; color=lift(a -> (:seagreen, a), addalpha))

    # --- result panel ---
    xl, yl = result_labels(expt)
    ax2 = Axis(left[2, 1]; xlabel=xl, ylabel=yl, title="Fit (active region)")
    gui[:ax_fit] = ax2
    hlines!(ax2, [0]; linewidth=0)
    errorbars!(ax2, state[:flat_errors]; whiskerwidth=8, color=state[:flat_error_colors])
    scatter!(ax2, state[:flat_points]; color=state[:flat_point_colors])
    lines!(ax2, state[:flat_fit]; color=state[:flat_fit_colors])
    text!(ax2, state[:seriestextpos]; text=state[:seriestexttxt],
          color=state[:seriestextcolor], fontsize=11, align=(:left, :center), offset=(6, 0))
    rowsize!(left, 1, Relative(0.55))
    on(_ -> autolimits!(ax2), state[:seriesdata])

    # --- controls ---
    r = 0
    r += 1
    right[r, 1] = Label(fig, "Spectrum:")
    sl = right[r, 2] = Slider(fig; range=1:state[:nplanes], width=170)
    connect!(state[:currentspectrum_idx], sl.value)

    if !singleregion
        r += 1
        right[r, 1] = Label(fig, "Region:")
        gui[:menu] = right[r, 2] = Menu(fig;
                                        options=lift(rs -> [reg.label for reg in rs],
                                                     state[:regions]), width=170)
        # Two-way sync between the dropdown and state[:active] (also driven by clicking or
        # hovering a region on the plot). Each side writes the other, so a reentrancy guard
        # is needed regardless of whether the underlying Observables suppress same-value
        # notifications - without it, a click bounces forever between the two handlers.
        syncingmenu = Ref(false)
        on(gui[:menu].i_selected) do i
            (isnothing(i) || syncingmenu[]) && return
            syncingmenu[] = true
            state[:active][] = i
            syncingmenu[] = false
        end
        on(state[:active]) do i
            (syncingmenu[] || !(1 ≤ i ≤ length(state[:regions][]))) && return
            syncingmenu[] = true
            gui[:menu].i_selected[] = i
            syncingmenu[] = false
        end
        # `on` only fires on future changes - set the initial displayed selection explicitly.
        gui[:menu].i_selected[] = state[:active][]
    end

    r += 1
    right[r, 1] = Label(fig, "Region width (ppm):")
    initw = isempty(state[:regions][]) ? 0.05 :
            state[:regions][][clamp(state[:active][], 1, length(state[:regions][]))].w
    tw = right[r, 2] = Textbox(fig; validator=Float64, width=170,
                               stored_string=string(round(initw; digits=3)))
    on(tw.stored_string) do s
        return set_region_width!(state, state[:active][], parse(Float64, s))
    end
    on(state[:active]) do i
        (1 ≤ i ≤ length(state[:regions][])) || return
        return tw.displayed_string[] = string(round(state[:regions][][i].w; digits=3))
    end

    r += 1
    gui[:togglefit] = Toggle(fig; active=true)
    right[r, 1] = gui[:togglefit]
    right[r, 2] = Label(fig, "Fitting")
    connect!(state[:isfitting], gui[:togglefit].active)

    r += 1
    right[r, 1] = Button(fig; label="(R)ename region")
    on(_ -> beginrename!(state), right[r, 1].clicks)
    right[r, 2] = Button(fig; label="(D)elete region")
    on(_ -> deleteregion!(state, state[:active][]), right[r, 2].clicks)

    r += 1
    right[r, 1:2] = Label(fig,
                          "Press A to add a region (tap for default width, drag to size it). Drag the shaded regions or the noise marker to reposition.";
                          word_wrap=true, tellwidth=false, halign=:left)

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

    setup_mouse!(fig, ax, state)
    setup_keyboard!(fig, ax, state)

    display(fig)
    autolimits!(ax2)
    while isopen(fig.scene)
        sleep(0.1)
    end
    GLMakie.closeall()
    return state
end

"""(Re)create one `vspan!` per region, coloured to highlight the active one. Called once
at startup and again whenever the number of regions changes (add/delete); simple
position/width/active-state changes are handled by the spans' own Observables and don't
need a rebuild."""
function rebuild_regionspans!(ax, state)
    gui = state[:gui]
    for sp in gui[:regionspans]
        delete!(ax, sp)
    end
    empty!(gui[:regionspans])
    for i in eachindex(state[:regions][])
        lo = lift(rs -> rs[i].c - rs[i].w / 2, state[:regions])
        hi = lift(rs -> rs[i].c + rs[i].w / 2, state[:regions])
        color = lift(a -> a == i ? ACTIVE_REGION_COLOR : INACTIVE_REGION_COLOR,
                     state[:active])
        push!(gui[:regionspans], vspan!(ax, lo, hi; color=color))
    end
end

"""
    beginrename!(state; initiator=:mouse)

Begin renaming the active region (shared by the button and the 'R' key). When triggered
`:keyboard`, the same keypress that opened rename mode also emits a unicode-input event
for its own character (e.g. 'r') a moment later — enter the intermediate `:renamingstart`
mode so `setup_keyboard!`'s unicode handler can swallow that spillover character rather
than prepending it to the new name (mirrors the 2D fitting GUI's rename handshake).
"""
function beginrename!(state; initiator=:mouse)
    state[:active][] > 0 || return
    state[:oldlabel][] = state[:activelabel][]
    state[:mode][] = initiator == :keyboard ? :renamingstart : :renaming
    return setactivelabel!(state, "‸")
end

"""Wire mouse dragging of the region/noise spans, and hover-highlighting of the region
under the cursor, on the spectrum axis."""
function setup_mouse!(fig, ax, state)
    gui = state[:gui]
    gui[:dragging] = :nothing
    on(events(ax).mousebutton) do ev
        ev.button == Mouse.left || return Consume(false)
        if ev.action == Mouse.press
            state[:mode][] == :normal || return Consume(false)
            gui[:dragging] = :nothing
            for (i, sp) in enumerate(gui[:regionspans])
                if mouseover(fig, sp)
                    gui[:dragging] = (:region, i)
                    state[:active][] = i
                    break
                end
            end
            if gui[:dragging] == :nothing && mouseover(fig, gui[:noisehit])
                gui[:dragging] = :noise
            end
            return Consume(gui[:dragging] != :nothing)
        elseif ev.action == Mouse.release
            gui[:dragging] = :nothing
            return Consume(false)
        end
    end
    on(events(fig).mouseposition; priority=2) do _
        if state[:mode][] == :addingdrag
            state[:addcurrent][] = mouseposition(ax)[1]
            return Consume(true)
        end
        state[:mode][] == :normal || return Consume(false)
        d = gui[:dragging]
        if d == :noise
            state[:noisec][] = mouseposition(ax)[1]
            return Consume(true)
        elseif d isa Tuple && d[1] == :region
            set_region_center!(state, d[2], mouseposition(ax)[1])
            return Consume(true)
        elseif d == :nothing
            # hover highlight: select whichever region the mouse is over (does not
            # affect dragging, which takes over as soon as the button is pressed)
            for (i, sp) in enumerate(gui[:regionspans])
                if mouseover(fig, sp)
                    state[:active][] = i
                    return Consume(false)
                end
            end
        end
        return Consume(false)
    end
end

"""Wire keyboard shortcuts: A to add a region (tap or drag), R to rename, D to delete,
Left/Right to step through spectra, matching the 2D fitting GUI's key bindings."""
function setup_keyboard!(fig, ax, state)
    on(events(fig).keyboardbutton; priority=2) do ev
        mode = state[:mode][]
        if mode == :normal && ev.action == Keyboard.press
            if ev.key == Keyboard.a
                state[:mode][] = :addingdrag
                x = mouseposition(ax)[1]
                state[:addstart][] = x
                state[:addcurrent][] = x
                return Consume(true)
            elseif ev.key == Keyboard.r
                beginrename!(state; initiator=:keyboard)
                return Consume(true)
            elseif ev.key == Keyboard.d
                deleteregion!(state, state[:active][])
                return Consume(true)
            elseif ev.key == Keyboard.left
                i = state[:currentspectrum_idx][]
                i > 1 && (state[:currentspectrum_idx][] = i - 1)
                return Consume(true)
            elseif ev.key == Keyboard.right
                i = state[:currentspectrum_idx][]
                i < state[:nplanes] && (state[:currentspectrum_idx][] = i + 1)
                return Consume(true)
            end
        elseif mode == :addingdrag
            if ev.action == Keyboard.release && ev.key == Keyboard.a
                lo, hi = minmax(state[:addstart][], state[:addcurrent][])
                if hi - lo < ADD_TAP_THRESHOLD
                    c = state[:addstart][]
                    addregion!(state, c - 0.025, c + 0.025)
                else
                    addregion!(state, lo, hi)
                end
                state[:mode][] = :renaming
                state[:oldlabel][] = state[:activelabel][]
                setactivelabel!(state, "‸")
                return Consume(true)
            end
        elseif mode == :renaming || mode == :renamingstart
            if ev.action == Keyboard.press
                label = state[:activelabel][]
                if ev.key == Keyboard.enter
                    setactivelabel!(state, strip(label[1:(end - 1)]))
                    state[:mode][] = :normal
                    return Consume(true)
                elseif ev.key == Keyboard.backspace
                    length(label) > 1 && setactivelabel!(state, label[1:(end - 2)] * "‸")
                    return Consume(true)
                elseif ev.key == Keyboard.escape
                    setactivelabel!(state, state[:oldlabel][])
                    state[:mode][] = :normal
                    return Consume(true)
                end
            end
            return Consume(false)
        end
        return Consume(false)
    end
    on(events(fig).unicode_input) do character
        if state[:mode][] == :renamingstart
            # the keypress that opened rename mode (the 'r' shortcut) also emits its own
            # character here a moment later - swallow it rather than prepending it
            state[:mode][] = :renaming
            character == 'r' && return Consume(true)
        elseif state[:mode][] != :renaming
            return Consume(false)
        end
        character in ('\r', '\n') && return Consume(false)  # handled as Keyboard.enter above
        label = state[:activelabel][]
        setactivelabel!(state, label[1:(end - 1)] * character * "‸")
        return Consume(true)
    end
end

"""
    pickregion(traces; peakppm=nothing, noiseppm=nothing, ppmwidth=0.05) -> (; peakppm, noiseppm, ppmwidth)
    pickregion(specs; kwargs...)

Interactively choose a single integration region and a noise position, returning the
integration triple. All traces are overlaid, so a region can be chosen that suits every
spectrum at once - which is what is wanted when one region must serve a whole set of
experiments (e.g. `Exchange1D`'s `integrate!`, which applies one triple to every
experiment in the problem).

This is the region-selection front-end on its own, without any fitting: use it where the
heavy fitting lives elsewhere.
"""
function pickregion(traces::AbstractVector{Trace}; peakppm=nothing, noiseppm=nothing,
                    ppmwidth=0.05)
    GLMakie.activate!(; focus_on_show=true,
                      title="NMRAnalysis.jl: select integration region")

    t1 = first(traces)
    peakppm = isnothing(peakppm) ? t1.δ[argmax(abs.(t1.y))] : peakppm
    noiseppm = isnothing(noiseppm) ? t1.δ[max(1, round(Int, 0.9 * length(t1.δ)))] : noiseppm

    pk = Observable(Float64(peakppm))
    nz = Observable(Float64(noiseppm))
    w = Observable(Float64(ppmwidth))

    fig = Figure(; size=(1000, 620))
    ax = Axis(fig[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title="Drag the peak region and noise marker, then press Accept")
    hlines!(ax, [0]; color=:grey)
    lines!(ax, _overlay_points([Point2f.(t.δ, t.y) for t in traces]); color=(:grey, 0.5),
           label="All spectra")
    peakspan = vspan!(ax, lift((p, ww) -> p - ww / 2, pk, w),
                      lift((p, ww) -> p + ww / 2, pk, w); color=ACTIVE_REGION_COLOR,
                      label="Peak")
    vlines!(ax, nz; color=:orchid, linewidth=2, label="Noise")
    noisehit = vspan!(ax, lift(c -> c - 0.01, nz), lift(c -> c + 0.01, nz); color=NOISE_COLOR)
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
            elseif mouseover(fig, noisehit)
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

"""Save the active-region fit figure and a text summary to the output folder."""
function save_results(state)
    dir = joinpath(pwd(), state[:outputdir][])
    isdir(dir) || mkpath(dir)
    expt = state[:expt]

    sd = state[:seriesdata][]
    xl, yl = result_labels(expt)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=xl, ylabel=yl)
    hlines!(ax, [0]; linewidth=0)
    for (i, s) in enumerate(sd)
        c = seriescolor(i)
        errorbars!(ax, s.errors; whiskerwidth=8, color=c)
        scatter!(ax, s.points; color=c, label=isempty(s.label) ? "Observed" : s.label)
        isempty(s.fitline) || lines!(ax, s.fitline; color=c)
    end
    length(sd) > 1 && axislegend(ax; position=:rt)
    save(joinpath(dir, "fit.pdf"), fig; backend=CairoMakie)

    open(joinpath(dir, "summary.txt"), "w") do f
        return println(f, state[:summary][])
    end
    @info "Saved results to $dir"
    return dir
end
