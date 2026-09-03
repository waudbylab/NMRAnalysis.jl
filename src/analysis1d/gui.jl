# Interactive single-window GUI for the 1D experiments. Mirrors the idioms of the R1ρ and
# GUI2D GUIs (GLMakie figure, draggable shaded regions, spectrum slider, live `lift`
# refit, keyboard add/rename/delete matching the 2D peak-picking interface). A second
# (exchange) window is deliberately out of scope here — none of relaxation / TRACT /
# nutation / diffusion / STD / kinetics needs it.

# A tap vs. a drag when adding a region with 'A' - shorter moves than this (ppm) are
# treated as a single click, giving the default-width region.
const ADD_TAP_THRESHOLD = 0.005

"""Scale the spectrum axis's y-view by `factor` - `factor=2` makes peaks look twice as
tall (divides both bounds by 2, i.e. halves the range), `factor=0.5` half as tall
(doubles the range). Scales `ymin`/`ymax` directly rather than re-centring on the
current view's midpoint, so it stays anchored to the data's own zero rather than
drifting off-centre over repeated presses. The axis's y-interaction is locked
(`yzoomlock`), so this is the controlled replacement for scroll/drag y-zoom - a
button/keyboard shortcut rather than the mouse."""
function scaleyaxis!(ax, factor)
    fl = ax.finallimits[]
    ylo, yhi = fl.origin[2], fl.origin[2] + fl.widths[2]
    return ylims!(ax, ylo / factor, yhi / factor)
end

"""Step the spectrum slider by `delta` planes (clamped to range), moving its visual
thumb via `set_close_to!` (not just `.value`, which the slider's thumb doesn't track -
see `setup_keyboard!`). Shared by the ←/→ buttons and the Left/Right key shortcuts."""
function stepslice!(sl, nplanes, delta)
    i = sl.value[] + delta
    return (1 ≤ i ≤ nplanes) && set_close_to!(sl, i)
end

"""
    gui!(expt::Experiment1D)

Launch the interactive analysis window for `expt`. The left column overlays all spectral
planes with draggable integration region(s) and a noise marker; the result panel below
shows the live fit for the active region. Returns the GUI state when the window closes.

# Keyboard shortcuts (mirroring the 2D fitting GUI)
- `A`: add a region. A quick tap adds a default-width region under the cursor; press,
  drag, and release to choose the width.
- `R`: rename the active region.
- `D`: delete the active region.
- Left/Right arrows: step through spectra.
- Up/Down arrows: scale the spectrum's y-axis (×2 / ÷2).
"""
function gui!(expt::Experiment1D)
    # the analysis type is already shown as a large bold label inside the window, so the
    # OS titlebar just carries the application identity rather than repeating it
    GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl")
    state = prepare_state(expt)
    state[:gui] = Dict{Symbol,Any}()
    gui = state[:gui]
    cols = Makie.wong_colors()
    singleregion = length(state[:regions][]) ≤ 1

    fig = Figure(; size=(1200, 800))
    left = fig[1, 1] = GridLayout()
    # Fixed, narrow width: `right`'s content is short labels/buttons, not something that
    # should compete with the spectra for space. Fixing it also makes every nested
    # GridLayout inside it (the button rows, the rename/delete pair) resolve against a
    # known, unambiguous width instead of an auto-determined one that shifts around as
    # content changes. `valign=:top` keeps its rows packed at the top, so any leftover
    # vertical space collects at the bottom instead of being spread between rows.
    right = fig[1, 2] = GridLayout(; valign=:top)
    colsize!(fig.layout, 2, Fixed(250))
    # The single outer row defaults to `Auto()`, which - since the axes report no
    # determinable height but the controls column's stack of labels/widgets does - was
    # being sized to the controls column's (shorter) natural height instead of the full
    # window, leaving blank space below. Force it to the whole available height.
    rowsize!(fig.layout, 1, Relative(1))

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
    # y-axis locked (matches the R1rho GUI): chemical-shift regions are picked/dragged
    # horizontally only, so zoom/pan/rectangle-select should never touch intensity. The
    # title identifies which plane is on display (e.g. "0.4 s delay (anti-TROSY)").
    spectitle = lift(i -> spectruminfo(expt, state[:planes].vars[i]), state[:currentspectrum_idx])
    ax = Axis(left[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title=spectitle, yzoomlock=true, ypanlock=true, xrectzoom=true, yrectzoom=false,
              xgridvisible=false, ygridvisible=false)
    gui[:ax_spec] = ax
    hlines!(ax, [0]; color=:grey)
    lines!(ax, state[:overlay]; color=(:grey, 0.35), label="All spectra")
    lines!(ax, state[:currentspectrum]; color=cols[1], label="Current spectrum")

    gui[:regionspans] = Any[]
    rebuild_regionspans!(ax, state)
    on(rs -> length(rs) != length(gui[:regionspans]) && rebuild_regionspans!(ax, state),
       state[:regions])

    # noise marker: a vline for visibility, plus a wider, mostly-transparent band as the
    # drag target (a true zero-width line is fiddly to grab with the mouse) - the larger
    # of 2% of the full spectral width or the widest region, so it's always at least as
    # easy to grab as the regions it estimates noise for, and reactive since regions can
    # be resized after the fact.
    basedefaultwidth = defaultregionwidth(first(state[:planes].traces).δ)
    noisehandlehw = lift(state[:regions]) do rs
        widest = isempty(rs) ? 0.0 : maximum(width(r) for r in rs)
        return max(basedefaultwidth, widest) / 2
    end
    vlines!(ax, state[:noisec]; color=:orchid, linewidth=2, label="Noise")
    gui[:noisehit] = vspan!(ax, lift((c, hw) -> c - hw, state[:noisec], noisehandlehw),
                            lift((c, hw) -> c + hw, state[:noisec], noisehandlehw);
                            color=NOISE_COLOR)
    axislegend(ax; position=:lt)

    # add-region drag preview
    addlo = lift((m, a, c) -> m == :addingdrag ? min(a, c) : 0.0, state[:mode],
                 state[:addstart], state[:addcurrent])
    addhi = lift((m, a, c) -> m == :addingdrag ? max(a, c) : 0.0, state[:mode],
                 state[:addstart], state[:addcurrent])
    addalpha = lift(m -> m == :addingdrag ? 0.25 : 0.0, state[:mode])
    # `xautolimits=false`: outside :addingdrag this span sits at the placeholder (0, 0),
    # and without opting out, that phantom point at ppm=0 gets pulled into the axis's
    # autolimits - forcing the view out to include zero even when the actual spectral
    # data never goes near it.
    vspan!(ax, addlo, addhi; color=lift(a -> (:seagreen, a), addalpha), xautolimits=false)

    # --- result panel ---
    # Built by the experiment's visualisation strategy; `gui[:ax_fit]` is set there.
    gui[:panelresult] = left[2, 1]
    resultpanel!(gui, state, expt)

    rowsize!(left, 1, Relative(0.55))

    # --- controls ---
    # `right` uses a single column throughout: every row is either one widget, or its
    # own small nested GridLayout for anything needing more than one widget side by
    # side. Sharing `right`'s own columns directly across differently-shaped rows (a
    # 7-item button row, a 2-item button pair, a label+textbox pair, ...) made Makie
    # stretch each row's cells to line up with whichever row needed the most columns,
    # which is what produced a very spread-out, misaligned layout. Now that `right`
    # itself is a fixed width (see above), each nested row-grid resolves against that
    # same known width instead of an ambiguous auto-determined one.
    r = 0
    r += 1
    right[r, 1] = Label(fig, "NMRAnalysis: $(windowtitle(expt))"; font=:bold, fontsize=16,
                        halign=:left, tellwidth=false)

    # zoom row: y-scale and reset only - slice navigation is its own row below, so each
    # row's few items comfortably fit the panel's fixed width
    r += 1
    zoomrow = right[r, 1] = GridLayout()
    btn_yup = zoomrow[1, 1] = Button(fig; label="×2 (↑)")
    btn_ydown = zoomrow[1, 2] = Button(fig; label="÷2 (↓)")
    on(_ -> scaleyaxis!(ax, 2), btn_yup.clicks)
    on(_ -> scaleyaxis!(ax, 0.5), btn_ydown.clicks)
    btn_resetzoom = zoomrow[1, 3] = Button(fig; label="reset zoom")
    # Clearing `ax.limits` (rather than calling `reset_limits!` directly) is what
    # actually gives a *fresh* full-data view: `ylims!`/`rebuild_regionspans!` pin
    # `ax.limits` to specific bounds as they run, and `reset_limits!` on its own just
    # re-applies whatever is currently pinned there rather than recomputing from data.
    on(_ -> (ax.limits[] = (nothing, nothing)), btn_resetzoom.clicks)

    # slice navigation row
    r += 1
    slicerow = right[r, 1] = GridLayout()
    btn_left = slicerow[1, 1] = Button(fig; label="←")
    btn_right = slicerow[1, 2] = Button(fig; label="→")
    sl = slicerow[1, 3] = Slider(fig; range=1:state[:nplanes], width=90)
    gui[:slider] = sl
    connect!(state[:currentspectrum_idx], sl.value)
    on(_ -> stepslice!(sl, state[:nplanes], -1), btn_left.clicks)
    on(_ -> stepslice!(sl, state[:nplanes], 1), btn_right.clicks)
    slicerow[1, 4] = Label(fig, lift(i -> "$i of $(state[:nplanes])", state[:currentspectrum_idx]))

    r += 1
    right[r, 1] = Label(fig,
                        "Press A to add a region (tap for default width, drag to size it).";
                        word_wrap=true, tellwidth=false, halign=:left)

    r += 1
    right[r, 1] = Label(fig,
                        "Drag the middle of a shaded region to move it, its edge to resize it, or the noise marker to reposition it.";
                        word_wrap=true, tellwidth=false, halign=:left)

    r += 1
    editrow = right[r, 1] = GridLayout()
    # explicit matching widths rather than relying on grid-stretch semantics to split
    # the row evenly - simpler and unambiguous
    btn_rename = editrow[1, 1] = Button(fig; label="(R)ename", width=115)
    on(btn_rename.clicks) do _
        return state[:mode][] == :normal && beginrename!(state)
    end
    btn_delete = editrow[1, 2] = Button(fig; label="(D)elete", width=115)
    on(btn_delete.clicks) do _
        return state[:mode][] == :normal && deleteregion!(state, state[:active][])
    end

    r += 1
    outputrow = right[r, 1] = GridLayout()
    # empty (not pre-filled with "out") so the placeholder text is visible - `state[:outputdir]`
    # already defaults to "out" independently (see `prepare_state`), so leaving this
    # untouched still saves there
    tout = outputrow[1, 1] = Textbox(fig; width=90, placeholder="out")
    on(tout.stored_string) do s
        return state[:outputdir][] = s
    end
    btn_save = outputrow[1, 2] = Button(fig; label="Save results")
    on(btn_save.clicks) do _
        return save_results(state)
    end

    r += 1
    fitrow = right[r, 1] = GridLayout()
    gui[:togglefit] = fitrow[1, 1] = Toggle(fig; active=true)
    fitrow[1, 2] = Label(fig, "Fitting")
    connect!(state[:isfitting], gui[:togglefit].active)
    fittingrow = r  # extra gap added below this row, once later rows exist to gap against

    if !singleregion
        r += 1
        regionrow = right[r, 1] = GridLayout()
        regionrow[1, 1] = Label(fig, "Region:")
        gui[:menu] = regionrow[1, 2] = Menu(fig;
                                            options=lift(rs -> [reg.label for reg in rs],
                                                         state[:regions]), width=140)
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
    right[r, 1] = Label(fig,
                        lift(lbl -> isempty(lbl) ? "No region selected" : "Region: $lbl",
                             state[:activelabel]);
                        font=:bold, halign=:left, tellwidth=false)

    r += 1
    right[r, 1] = Label(fig, state[:resultsheader]; tellwidth=false, halign=:left,
                        justification=:left)

    r += 1
    right[r, 1] = Label(fig, state[:secondaryresult]; tellwidth=false, halign=:left,
                        justification=:left, word_wrap=true)

    # now that every later row exists, add breathing room below the fitting toggle
    rowgap!(right, fittingrow, Fixed(20))

    setup_mouse!(fig, ax, state)
    setup_keyboard!(fig, ax, state)

    display(fig)
    autolimits!(gui[:ax_fit])
    while isopen(fig.scene)
        sleep(0.1)
    end
    GLMakie.closeall()
    return state
end

"""(Re)create one `vspan!` per region, coloured to highlight the active one. Called once
at startup and again whenever the number of regions changes (add/delete); simple
position/width/active-state changes are handled by the spans' own Observables and don't
need a rebuild.

Each span's position/colour Observables are plain `Observable`s kept in sync by an
explicit `on` listener (rather than `lift`ing straight off `state[:regions]`/
`state[:active]`, index `i` baked in) so the listener can be `off`'d here before the next
rebuild. Without that, an old listener for index `i` stays registered on `state[:regions]`
after its span is deleted, and fires (indexing out of bounds) the next time a region is
removed and the list is shorter than `i`."""
function rebuild_regionspans!(ax, state)
    gui = state[:gui]
    # Capture the current view and reapply it once the spans below are rebuilt. Deleting
    # and re-adding plot objects makes Makie fall back to full-data autolimits on the
    # next relimit, discarding whatever the user had zoomed/panned to. Only on a genuine
    # rebuild (spans already exist): before the first render, `finallimits` is still its
    # placeholder default rather than the real data range, so capturing it on the
    # initial build would lock the axis to that tiny placeholder instead of letting the
    # normal first-time autolimits happen.
    #
    # Restored via `targetlimits` (the transient "what's currently shown" state), not
    # `xlims!`/`ax.limits` (a *permanent* override): an earlier version used `xlims!`
    # here, which pins `ax.limits` until something explicitly clears it - and then the
    # "reset zoom" button and the ×2/÷2 y-scale buttons (which read `ax.limits` back via
    # `reset_limits!`/`ylims!` to decide what to leave unchanged) would silently re-apply
    # that stale pin, snapping the x-view back to whatever it was at the last add/delete
    # instead of leaving it alone or truly resetting it.
    isrebuild = !isempty(gui[:regionspans])
    preservedview = ax.finallimits[]

    for (sp, obsfuncs) in gui[:regionspans]
        delete!(ax, sp)
        foreach(off, obsfuncs)
    end
    empty!(gui[:regionspans])
    for i in eachindex(state[:regions][])
        lo = Observable(0.0)
        hi = Observable(0.0)
        color = Observable(INACTIVE_REGION_COLOR)
        of1 = on(state[:regions]; update=true) do rs
            i <= length(rs) || return
            r = rs[i]
            lo[] = r.lo
            hi[] = r.hi
        end
        of2 = on(state[:active]; update=true) do a
            color[] = a == i ? ACTIVE_REGION_COLOR : INACTIVE_REGION_COLOR
        end
        sp = vspan!(ax, lo, hi; color=color)
        push!(gui[:regionspans], (sp, (of1, of2)))
    end
    isrebuild && (ax.targetlimits[] = preservedview)
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
            for (i, (sp, _)) in enumerate(gui[:regionspans])
                if mouseover(fig, sp)
                    # dragging the outer 25% of the region drags that edge (moving both
                    # centre and width); the middle 50% drags the whole region instead
                    r = state[:regions][][i]
                    lo, hi = r.lo, r.hi
                    w = width(r)
                    frac = w > 0 ? (mouseposition(ax)[1] - lo) / w : 0.5
                    gui[:dragging] = if frac < 0.25
                        (:regionedge, i, hi)
                    elseif frac > 0.75
                        (:regionedge, i, lo)
                    else
                        (:region, i)
                    end
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
        elseif d isa Tuple && d[1] == :regionedge
            set_region_edge!(state, d[2], d[3], mouseposition(ax)[1])
            return Consume(true)
        elseif d == :nothing
            # hover highlight: select whichever region the mouse is over (does not
            # affect dragging, which takes over as soon as the button is pressed)
            for (i, (sp, _)) in enumerate(gui[:regionspans])
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
Left/Right to step through spectra, Up/Down to scale the spectrum's y-axis, matching the
2D fitting GUI's key bindings."""
function setup_keyboard!(fig, ax, state)
    defaultwidth = defaultregionwidth(first(state[:planes].traces).δ)
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
                stepslice!(state[:gui][:slider], state[:nplanes], -1)
                return Consume(true)
            elseif ev.key == Keyboard.right
                stepslice!(state[:gui][:slider], state[:nplanes], 1)
                return Consume(true)
            elseif ev.key == Keyboard.up
                scaleyaxis!(ax, 2)
                return Consume(true)
            elseif ev.key == Keyboard.down
                scaleyaxis!(ax, 0.5)
                return Consume(true)
            end
        elseif mode == :addingdrag
            if ev.action == Keyboard.release && ev.key == Keyboard.a
                lo, hi = minmax(state[:addstart][], state[:addcurrent][])
                if hi - lo < ADD_TAP_THRESHOLD
                    c = state[:addstart][]
                    addregion!(state, c - defaultwidth / 2, c + defaultwidth / 2)
                else
                    addregion!(state, lo, hi)
                end
                state[:mode][] = :normal
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
    pickregion(traces; peakppm=nothing, noiseppm=nothing, ppmwidth=defaultregionwidth(...)) -> (; peakppm, noiseppm, ppmwidth)
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
                    ppmwidth=defaultregionwidth(first(traces).δ))
    GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl")

    t1 = first(traces)
    peakppm = isnothing(peakppm) ? t1.δ[argmax(abs.(t1.y))] : peakppm
    noiseppm = isnothing(noiseppm) ? t1.δ[max(1, round(Int, 0.9 * length(t1.δ)))] : noiseppm

    pk = Observable(Float64(peakppm))
    nz = Observable(Float64(noiseppm))
    w = Observable(Float64(ppmwidth))

    fig = Figure(; size=(1000, 620))
    # single column: force it to the full window width, otherwise it shrinks to the
    # (narrower, content-determined) control row's width instead of filling the window
    colsize!(fig.layout, 1, Relative(1))
    ax = Axis(fig[1, 1]; xreversed=true, xlabel="Chemical shift (ppm)", ylabel="Intensity",
              title="Drag the peak region and noise marker, then press Accept",
              yzoomlock=true, ypanlock=true, xrectzoom=true, yrectzoom=false,
              xgridvisible=false, ygridvisible=false)
    hlines!(ax, [0]; color=:grey)
    # plot at most the first 10 traces - with e.g. a whole titration or pseudo-3D series,
    # overlaying every one of them is slow to render and no more informative than a
    # representative handful
    plotted = traces[1:min(10, length(traces))]
    lines!(ax, _overlay_points([Point2f.(t.δ, t.y) for t in plotted]); color=(:grey, 0.5),
           label="All spectra")
    peakspan = vspan!(ax, lift((p, ww) -> p - ww / 2, pk, w),
                      lift((p, ww) -> p + ww / 2, pk, w); color=ACTIVE_REGION_COLOR,
                      label="Peak")
    vlines!(ax, nz; color=:orchid, linewidth=2, label="Noise")
    noisehandlehw = defaultregionwidth(t1.δ) / 2
    noisehit = vspan!(ax, lift(c -> c - noisehandlehw, nz), lift(c -> c + noisehandlehw, nz);
                      color=NOISE_COLOR)
    axislegend(ax; position=:lt)

    # restrict the view to the width common to every input trace, rather than the union -
    # with several experiments of differing sweep width, the wider ones would otherwise
    # pad the view with a region where some inputs have no data at all
    lo = maximum(minimum(t.δ) for t in traces)
    hi = minimum(maximum(t.δ) for t in traces)
    # `xlims!(ax, a, b)` sets `ax.xreversed[] = a > b`, so the larger value has to come
    # first here too, to keep this axis's chemical-shift convention.
    lo < hi && xlims!(ax, hi, lo)

    ctrl = fig[2, 1] = GridLayout()
    ctrl[1, 1] = Label(fig, "Integration width (ppm):")
    tb = ctrl[1, 2] = Textbox(fig; stored_string=string(round(w[]; digits=3)),
                              validator=Float64, width=100)
    on(s -> w[] = parse(Float64, s), tb.stored_string)
    on(neww -> tb.displayed_string[] = string(round(neww; digits=3)), w)
    accept = ctrl[1, 3] = Button(fig; label="Accept")
    done = Ref(false)
    on(_ -> done[] = true, accept.clicks)

    dragging = Ref{Any}(:nothing)
    on(events(ax).mousebutton) do ev
        ev.button == Mouse.left || return Consume(false)
        if ev.action == Mouse.press
            dragging[] = if mouseover(fig, peakspan)
                # dragging the outer 25% of the span drags that edge (moving both
                # centre and width, mirroring the main GUI's region-resize); the
                # middle 50% drags the whole span instead
                lo, hi = pk[] - w[] / 2, pk[] + w[] / 2
                frac = w[] > 0 ? (mouseposition(ax)[1] - lo) / w[] : 0.5
                if frac < 0.25
                    (:edge, hi)
                elseif frac > 0.75
                    (:edge, lo)
                else
                    :peak
                end
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
        d = dragging[]
        if d == :peak
            pk[] = mouseposition(ax)[1]
            return Consume(true)
        elseif d == :noise
            nz[] = mouseposition(ax)[1]
            return Consume(true)
        elseif d isa Tuple && d[1] == :edge
            fixed = d[2]
            lo, hi = minmax(fixed, mouseposition(ax)[1])
            pk[] = (lo + hi) / 2
            w[] = max(hi - lo, 1e-6)
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

"""Save the fit figure(s) and a text summary, covering every region (not just the active
one shown live in the GUI panels), to the output folder: an overlay of every region in
`fit.pdf`, plus one `fit_<region>.pdf` per region on its own."""
function save_results(state)
    dir = joinpath(pwd(), state[:outputdir][])
    isdir(dir) || mkpath(dir)
    expt = state[:expt]
    result = state[:result][]
    labels = [r.label for r in state[:regions][]]
    xl, yl = resultlabels(expt)

    # no gridlines (matching the live GUI panels), but keep a visible zero line
    newfitaxis(fig) = Axis(fig[1, 1]; xlabel=xl, ylabel=yl, xgridvisible=false,
                           ygridvisible=false)

    # overlay of every region
    fig = Figure()
    ax = newfitaxis(fig)
    hlines!(ax, [0]; color=:grey)
    i = 0
    for label in labels
        i += plotresult!(ax, expt, result, label, i)
    end
    i > 1 && axislegend(ax; position=:rt)
    save(joinpath(dir, "fit.pdf"), fig; backend=CairoMakie)

    # one plot per region
    for label in labels
        fig1 = Figure()
        ax1 = newfitaxis(fig1)
        hlines!(ax1, [0]; color=:grey)
        n1 = plotresult!(ax1, expt, result, label, 0)
        n1 > 1 && axislegend(ax1; position=:rt)
        save(joinpath(dir, "fit_$label.pdf"), fig1; backend=CairoMakie)
    end

    open(joinpath(dir, "summary.txt"), "w") do f
        println(f, "Integration regions (ppm):")
        for r in state[:regions][]
            println(f,
                    "  $(r.label): $(round(r.lo; digits=4)) to $(round(r.hi; digits=4))")
        end
        println(f, "Noise position (ppm): $(round(state[:noisec][]; digits=4))")
        println(f)
        for label in labels
            print(f, summary_text(expt, result, label))
        end
    end
    @info "Saved results to $dir"
    return dir
end
