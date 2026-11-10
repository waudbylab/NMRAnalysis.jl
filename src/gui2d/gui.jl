"""
    gui!(expt; peaklist=nothing, call=nothing) -> expt

Open the analysis window for `expt`, loading the peaks in `peaklist` if given (a peak list,
or a folder holding `peaklist.csv`), and wait until it is closed, returning `expt` with its fitted peaks (see
[`results`](@ref) and [`planeresults`](@ref)). `call` is the [`AnalysisCall`](@ref) the
entry point recorded, written to `summary.txt` so the analysis can be repeated.
"""
function gui!(expt::Experiment; peaklist=nothing, call=nothing)
    GLMakie.activate!(; title=windowname("2D " * analysisname(expt)), focus_on_show=true)
    # No gridlines anywhere (interactive or exported plots use the same Makie theme); axes
    # still show a visible zero line where one is drawn explicitly.
    set_theme!(Theme(;
                     Axis=(xgridvisible=false, ygridvisible=false,
                           xminorgridvisible=false, yminorgridvisible=false)))

    state = expt.state[]
    state[:call] = Observable{Union{AnalysisCall,Nothing}}(call)

    g = Dict{Symbol,Any}() # GUI state
    state[:gui] = g

    g[:fig] = Figure(; size=(1200, 800))
    g[:paneltop] = g[:fig][1, 1:2] = GridLayout()
    g[:panelcontour] = g[:fig][2:3, 1] = GridLayout()
    g[:panelinfo] = g[:fig][2, 2] = GridLayout()
    g[:panelpeakplot] = g[:fig][3, 2] = GridLayout()
    g[:status] = Label(g[:fig][4, 1:2], statustext(expt); tellwidth=false, halign=:center)

    rowsize!(g[:fig].layout, 1, Auto(true))
    # Row 2 (info panel: buttons, hint, radius sliders) is sized to its own content so it
    # never gets squeezed; row 3 (peak plot) takes whatever height is left.
    rowsize!(g[:fig].layout, 2, Auto(true))
    rowsize!(g[:fig].layout, 3, Auto(false, 1))
    # Row 4, the status line, has a fixed height so that text coming and going doesn't
    # shift the panels above it.
    rowsize!(g[:fig].layout, 4, Fixed(24))
    colsize!(g[:fig].layout, 1, Auto(false, 1.6))
    colsize!(g[:fig].layout, 2, Auto(false, 1))

    # top panel
    g[:cmdcontourup] = Button(g[:paneltop][1, 1]; label="contour ↑")
    g[:cmdcontourdown] = Button(g[:paneltop][1, 2]; label="contour ↓")
    g[:cmdresetzoom] = Button(g[:paneltop][1, 3]; label="reset zoom")
    g[:sliderslice] = Slider(g[:paneltop][1, 4]; range=1:nslices(expt))
    g[:cmdsliceleft] = Button(g[:paneltop][1, 5]; label="←")
    g[:cmdsliceright] = Button(g[:paneltop][1, 6]; label="→")
    g[:slicelabel] = Label(g[:paneltop][1, 7], state[:current_slice_label])

    # Moving-peak experiments get a "Show all" planes toggle just left of the Fitting toggle.
    # The toggle's plots and handler are wired up later in add_moving_overlays!.
    col = 8
    if !hasfixedpositions(expt)
        Label(g[:paneltop][1, col], "(S)how all")
        g[:toggleother] = Toggle(g[:paneltop][1, col + 1]; active=true)
        col += 2
    end
    Label(g[:paneltop][1, col], "Fitting")
    g[:togglefit] = Toggle(g[:paneltop][1, col + 1]; active=true)
    g[:cmdsummary] = Button(g[:paneltop][1, col + 2]; label="Summary plot")
    g[:cmdquit] = Button(g[:paneltop][1, col + 3]; label="(Q)uit")

    # create contour plot
    g[:basecontour] = Observable(10.0)
    g[:contourscale] = Observable(1.7)
    g[:logscale] = lift(r -> [r .^ (0:10); -1 * (r .^ (0:10))], g[:contourscale])
    g[:contourlevels] = lift((c0, arr) -> c0 * arr, g[:basecontour], g[:logscale])

    g[:axcontour] = Axis(g[:panelcontour][1, 1];
                         xlabel="$(label(expt.specdata.nmrdata[1],F1Dim)) chemical shift (ppm)",
                         ylabel="$(label(expt.specdata.nmrdata[1],F2Dim)) chemical shift (ppm)",
                         xreversed=true, yreversed=true)
    g[:pltmask] = heatmap!(g[:axcontour],
                           state[:current_mask_x],
                           state[:current_mask_y],
                           state[:current_mask_z];
                           colormap=[:white, :lightgoldenrod1],
                           colorrange=(0, 1))
    g[:pltcontour] = contour!(g[:axcontour],
                              state[:current_spec_x],
                              state[:current_spec_y],
                              state[:current_spec_z];
                              levels=g[:contourlevels],
                              signcolours(:grey50, :lightblue)...)
    g[:pltfit] = contour!(g[:axcontour],
                          state[:current_fit_x],
                          state[:current_fit_y],
                          state[:current_fit_z];
                          levels=g[:contourlevels],
                          signcolours(:orangered, :dodgerblue)...)

    # # create 3D plot
    # g[:ax3d] = Axis3(g[:panel3d][1,1], xlabel="δX / ppm", ylabel="δy / ppm", zlabel="Intensity",
    #     xreversed=true, yreversed=true)
    # g[:pltwireobs] = wireframe!(g[:ax3d],
    #     state[:current_spec_x],
    #     state[:current_spec_y],
    #     state[:current_spec_z],
    #     color=:grey80)
    # g[:pltwirefit] = wireframe!(g[:ax3d],
    #     state[:current_fit_x],
    #     state[:current_fit_y],
    #     state[:current_fit_z],
    #     color=:orangered)

    # add the peak positions
    g[:pltpeaks] = scatter!(g[:axcontour], state[:positions]; markersize=10, marker=:x,
                            color=state[:peakcolours])
    g[:pltlabels] = text!(g[:axcontour], state[:positions]; text=state[:labels],
                          fontsize=14,
                          font=:bold,
                          offset=(8, 0),
                          align=(:left, :center),
                          color=:black)
    # The drag handle. The moused-over/selected peak is enlarged (see state.jl) for grab
    # feedback, and the handle is kept at the top z-level so it stays pickable above the
    # trajectory/context overlays.
    g[:pltinitialpeaks] = scatter!(g[:axcontour], state[:initialpositions];
                                   markersize=state[:initialpeaksizes],
                                   color=state[:peakcolours])
    translate!(g[:pltinitialpeaks], 0, 0, 10)

    # moving-peak overlays (per-plane position trajectories + faint context planes);
    # a no-op for fixed-position experiments
    add_moving_overlays!(g, state, expt)

    # peak info
    Label(g[:panelinfo][1, 1:2], "Working directory:\n$(shortpath(pwd()))"; word_wrap=true,
          tellwidth=false, halign=:left)
    g[:cmdload] = Button(g[:panelinfo][2, 1]; label="Load peak list")
    # Output folder typed rather than chosen from a native dialog: the dialog cannot create
    # a folder on every platform, and a typed name is what makes a save repeatable.
    outputrow = g[:panelinfo][2, 2] = GridLayout()
    g[:toutput] = Textbox(g[:fig]; width=90, stored_string=expt.state[][:outputdir][])
    outputrow[1, 1] = g[:toutput]
    # An empty box means the default rather than the working directory itself.
    on(g[:toutput].stored_string) do s
        return expt.state[][:outputdir][] = isempty(strip(s)) ? "out" : strip(s)
    end
    commitondefocus!(g[:toutput])
    g[:cmdsave] = Button(g[:fig]; label="Save")
    outputrow[1, 2] = g[:cmdsave]
    Label(g[:panelinfo][3, 1:2], addpeakhint(expt); fontsize=14, halign=:left,
          justification=:left, tellwidth=false)
    g[:cmdrename] = Button(g[:panelinfo][4, 1]; label="(R)ename peak")
    g[:cmddelete] = Button(g[:panelinfo][4, 2]; label="(D)elete peak")
    # A SliderGrid's slider column only expands into space its container already has, so
    # left to auto-size it collapses to a token width; giving it an explicit width is what
    # makes the track usable.
    g[:sgradiix] = SliderGrid(g[:panelinfo][5, 1],
                              (label="X radius", range=0.02:0.005:0.1, format="{:.3f} ppm",
                               startvalue=expt.xradius[]); width=190)
    g[:sgradiiy] = SliderGrid(g[:panelinfo][5, 2],
                              (label="Y radius", range=0.1:0.02:0.8, format="{:.2f} ppm",
                               startvalue=expt.yradius[]); width=190)
    g[:sliderxradius] = g[:sgradiix].sliders[1].value
    g[:slideryradius] = g[:sgradiiy].sliders[1].value

    g[:infotext] = Label(g[:panelinfo][6, 1:2], state[:current_peak_info])

    # peak plot panel
    makepeakplot!(g, state, expt)

    @debug "Adding handlers"
    addhanders!(g, state, expt)
    set_close_to!(g[:sliderslice], initialslice(expt))
    isnothing(peaklist) || loadpeaks!(expt, peaklist)

    capturewindow(g[:fig]; ready=() -> settle!(expt)) && return expt
    display(g[:fig])
    while isopen(g[:fig].scene)
        sleep(0.1)
    end
    return expt
end

"""
For a screenshot: wait for the fit in progress, set the contours, then select and zoom in
on the peak it asks for.
"""
function settle!(expt)
    state = expt.state[]
    g = state[:gui][]
    shot = SCREENSHOT[]
    while (task = state[:fit_task][]; !isnothing(task) && !istaskdone(task))
        wait(task)
    end
    g[:basecontour][] *= g[:contourscale][]^shot.contour
    peaks = expt.peaks[]
    isempty(peaks) && return nothing
    i = isnothing(shot.peak) ? 1 : findfirst(p -> p.label[] == shot.peak, peaks)
    isnothing(i) && throw(ArgumentError("no peak labelled $(shot.peak)"))
    state[:current_peak_idx][] = i
    if shot.zoom != 1
        x, y = state[:positions][][i]
        xlo, xhi = extrema(state[:current_spec_x][])
        ylo, yhi = extrema(state[:current_spec_y][])
        hx, hy = (xhi - xlo) / 2shot.zoom, (yhi - ylo) / 2shot.zoom
        # limits! sets each axis reversed when its first value is the larger, as shifts are
        limits!(g[:axcontour], x + hx, x - hx, y + hy, y - hy)
    end
    return nothing
end

"""
    statustext(expt) -> Observable{String}

The status line: what the window is doing or waiting for, the outcome of the last command,
and what is left to fit, with the key that fits it.
"""
function statustext(expt::Experiment)
    state = expt.state[]
    return lift(state[:mode], state[:fitprogress], state[:fitround], state[:message],
                state[:lineheld], expt.peaks,
                expt.isfitting) do mode, progress, r, message, held, peaks, fitting
        mode == :fitting && return fittingtext(progress, r)
        mode == :line && return linehint(held)
        haskey(MODEHINTS, mode) && return MODEHINTS[mode]
        return join(filter(!isempty, [message; pendingtext(peaks, fitting)]), "   ·   ")
    end
end

function fittingtext((done, total), round)
    progress = "fitting cluster $(min(done + 1, total)) of $total"
    round > 0 && return "Continuing to convergence, round $round: $progress (Esc to stop)"
    return uppercasefirst(progress) * " (Esc to stop)"
end

"One entry for each kind of peak whose fit is not finished, saying what finishes it."
function pendingtext(peaks, fitting)
    fitted = filter(peak -> !peak.touched[], peaks)
    waiting = length(peaks) - length(fitted)
    stopped = count(iscontinuable, fitted)
    atlimit = count(peak -> fitstatus(peak) == :bound, fitted)
    lines = String[]
    waiting > 0 && push!(lines,
                         counted(waiting, "peak") * " awaiting a fit: " *
                         (fitting ? "press (C) to fit" : "switch Fitting on to fit"))
    stopped > 0 && push!(lines,
                         counted(stopped, "peak") *
                         " stopped short: (C) continue, ⇧C to converge")
    atlimit > 0 && push!(lines,
                         counted(atlimit, "peak") *
                         " at a limit: move, or widen radii (⇧←↑↓→)")
    return lines
end

counted(n, noun) = n == 1 ? "1 $noun" : "$n $(noun)s"

# The keys each mode is waiting for. :moving needs none: the drag is self-explanatory.
const RENAMEHINT = "Renaming: type the new label   (Enter) finish   (Esc) cancel"
const MODEHINTS = Dict(:saving => "Saving results…",
                       :renaming => RENAMEHINT,
                       :renamingstart => RENAMEHINT,
                       :adding =>
                           "Adding a peak plane by plane:   (A) mark this plane   " *
                           "(Space) fill the rest   (Esc) cancel")

# Marking a line, by dragging with (L) held, or by two presses or a press and a click
function linehint(held)
    held && return "Dragging a line: release (L) at its other end   (Esc) cancel"
    return "Drawing a line: click its other end, or press (L) there   (Esc) cancel"
end

"""
    notice!(expt, text)

Show `text` in the status line until the peaks are next changed or fitted.
"""
notice!(expt::Experiment, text) = (expt.state[][:message][] = text)

"Say in the status line that the newest peak was added, and `how`."
function noticeadded!(expt::Experiment, how="")
    return notice!(expt, "Added peak $(expt.peaks[][end].label[])$how")
end

"Delete the selected peak, if there is one, and say which in the status line."
function deleteselected!(expt::Experiment)
    state = expt.state[]
    idx = state[:current_peak_idx][]
    idx > 0 || return nothing
    label = expt.peaks[][idx].label[]
    state[:current_peak_idx][] = 0
    deletepeak!(expt, idx)
    return notice!(expt, "Deleted peak $label")
end

# Hook for experiment-specific contour-panel overlays; specialised for moving-peak
# experiments in expt-moving.jl. No-op for fixed-position experiments.
add_moving_overlays!(g, state, ::Experiment) = nothing

# Shortcuts shown above the peak-info panel, one per line; moving-peak experiments also
# advertise (T)rack and (L)ine.
addpeakhint(::Experiment) = "(A) add a peak under the cursor\n" * RADIUSHINT

const RADIUSHINT = "⇧←↑↓→ set the selected peak's own radii"

function addhanders!(g, state, expt::Experiment)
    g[:fig].scene.backgroundcolor = lift(state[:mode]) do mode
        if mode == :fitting
            RGBAf(1.0, 0.63, 0.48, 1.0)      # :salmon
        elseif mode == :renaming || mode == :renamingstart
            RGBAf(0.75, 0.94, 1.0, 1.0)     # :lightblue
        elseif mode == :moving
            RGBAf(0.6, 0.98, 0.6, 1.0)      # :palegreen
        elseif mode == :adding || mode == :line
            RGBAf(1.0, 0.95, 0.7, 1.0)      # light yellow
        elseif mode == :saving
            RGBAf(0.8, 0.8, 0.8, 1.0)       # :grey80, busy signal while writing results
        else
            RGBAf(1.0, 1.0, 1.0, 1.0)       # :white
        end
    end

    # link 2D and 3D contour plot limits
    if haskey(g, :ax3d)
        on(g[:axcontour].xaxis.attributes.limits) do xl
            return xlims!(g[:ax3d], xl)
        end
        on(g[:axcontour].yaxis.attributes.limits) do yl
            return ylims!(g[:ax3d], yl)
        end
    end

    # contour levels
    on(g[:cmdcontourup].clicks) do _
        return g[:basecontour][] *= g[:contourscale][]
    end
    on(g[:cmdcontourdown].clicks) do _
        return g[:basecontour][] /= g[:contourscale][]
    end

    # zoom
    on(g[:cmdresetzoom].clicks) do _
        return reset_limits!(g[:axcontour])
    end

    # slices
    connect!(state[:current_slice], g[:sliderslice].value)
    on(g[:cmdsliceleft].clicks) do _
        i = state[:current_slice][]
        if i > 1
            i -= 1
            set_close_to!(g[:sliderslice], i)
        end
    end
    on(g[:cmdsliceright].clicks) do _
        i = state[:current_slice][]
        if i < nslices(expt)
            i += 1
            set_close_to!(g[:sliderslice], i)
        end
    end

    # fitting active
    connect!(expt.isfitting, g[:togglefit].active)
    on(x -> g[:pltfit].visible = x, g[:togglefit].active)

    # summary plot button: grey out when no peaks
    function _update_summary_button(peaks)
        if isempty(peaks)
            g[:cmdsummary].buttoncolor[] = RGBAf(0.85, 0.85, 0.85, 1.0)
            g[:cmdsummary].labelcolor[] = RGBAf(0.6, 0.6, 0.6, 1.0)
        else
            g[:cmdsummary].buttoncolor[] = RGBAf(0.94, 0.94, 0.94, 1.0)
            g[:cmdsummary].labelcolor[] = RGBAf(0.0, 0.0, 0.0, 1.0)
        end
    end
    _update_summary_button(expt.peaks[])  # set initial state
    on(_update_summary_button, expt.peaks)

    on(g[:cmdsummary].clicks) do _
        isempty(expt.peaks[]) && return
        @async display(GLMakie.Screen(), summaryplot(expt))
    end

    # load the peak list saved in the output folder
    on(g[:cmdload].clicks) do _
        folder = joinpath(pwd(), expt.state[][:outputdir][])
        try
            loadpeaks!(expt, folder)
        catch e
            @error "Error loading the peak list in $folder" exception = (e,
                                                                         catch_backtrace())
            notice!(expt, "Could not load the peak list: peaks unchanged, see the REPL")
        end
    end

    # save peak list
    on(g[:cmdsave].clicks) do _
        return saveresults!(expt)
    end

    # quit
    on(g[:cmdquit].clicks) do _
        return quit!()
    end

    # delete peak
    on(_ -> deleteselected!(expt), g[:cmddelete].clicks)

    # rename peak
    on(g[:cmdrename].clicks) do _
        if state[:current_peak_idx][] > 0
            renamepeak!(expt, state, :mouse)
        end
    end

    # radius sliders
    connect!(expt.xradius, g[:sliderxradius])
    connect!(expt.yradius, g[:slideryradius])

    # peak hover
    onpick(g[:axcontour], g[:pltinitialpeaks]) do _, idx
        # if state[:current_peak_idx][] != idx && state[:mode][] == :normal
        # GPU picking can lag a frame behind a peak deletion, so idx may momentarily point
        # past the end of the (now shorter) peak list - ignore a stale pick rather than crash.
        if state[:mode][] == :normal && idx in eachindex(expt.peaks[])
            @debug "Setting current peak to $idx"
            state[:current_peak_idx][] = idx
            if haskey(g, :axpeakplot)
                autolimits!(g[:axpeakplot])
            end
            if haskey(g, :axpeakplotX)
                autolimits!(g[:axpeakplotX])
            end
            if haskey(g, :axpeakplotY)
                autolimits!(g[:axpeakplotY])
            end
        end
    end
    # mouse handling
    on(events(g[:axcontour]).mousebutton; priority=2) do event
        return process_mousebutton(expt, state, event)
    end
    on(events(g[:axcontour]).mouseposition; priority=2) do mousepos
        return process_mouseposition(expt, state, mousepos)
    end

    # keyboard
    on(events(g[:axcontour]).keyboardbutton; priority=2) do event
        return process_keyboardbutton(expt, state, event)
    end
    on(events(g[:fig]).unicode_input) do character
        return process_unicode_input(expt, state, character)
    end

    # window cleanup on close. window_open fires `false` when the window is closed; without
    # explicit cleanup GLMakie can leave an invisible window frame behind (notably on macOS),
    # which otherwise requires a manual GLMakie.closeall(). We also cancel any in-flight fit so a
    # background task isn't left writing into a closed window's state.
    on(events(g[:fig].scene).window_open) do isopen
        isopen && return
        @debug "Window closed - cleaning up"
        state[:fit_generation][] += 1
        # closeall() clears the ghost frame. The app shows one GUI at a time, so closing all
        # GLMakie windows here is acceptable; switch to a targeted screen close if multiple
        # simultaneous windows are ever needed.
        #
        # This observer fires from inside the renderloop's poll cycle, so calling closeall()
        # synchronously tears down the GL context while the renderloop is still running it
        # ("Context is not alive anymore!"). Defer it to a separate task (as a manual REPL
        # closeall() would be) and let the closing window's own teardown settle first.
        @async begin
            sleep(0.2)
            GLMakie.closeall()
        end
        return
    end
end

function renamepeak!(expt, state, initiator)
    @debug "Renaming peak"
    if initiator == :keyboard
        state[:mode][] = :renamingstart
    else
        state[:mode][] = :renaming
    end
    state[:oldlabel][] = state[:current_peak][].label[]
    state[:current_peak][].label[] = "‸"
    return notify(expt.peaks)
end

"Close the GUI window, deferred so the click/keypress event finishes first."
function quit!()
    return @async begin
        sleep(0.05)
        GLMakie.closeall()
    end
end

# Colour contours by sign of level. Makie drops levels outside the data range,
# so a per-level colour vector no longer matches the number of levels.
function signcolours(pos, neg)
    return (colormap=[neg, pos], colorrange=(-1.0f-12, 1.0f-12), lowclip=neg, highclip=pos)
end
