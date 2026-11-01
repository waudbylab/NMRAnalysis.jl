# The interactive window.

const HELP = ("←→ / 0-9 select   r rename   d delete line   g group   u ungroup   " *
              "drag: trajectory   click: pick   drag marker: move   ⇧-drag: erase   " *
              "right-click: delete pick",
              "m snap on/off (⇧ inverts)   ↑↓ scale   [ ] radius   a align   ⌫ undo   " *
              "s save figure   q quit   scroll: zoom   right-drag: pan   ctrl-click: reset view")

"Gain step for the up and down keys."
const STEP = 2^(1 / 4)

shifted(scene) = ispressed(scene, Keyboard.left_shift | Keyboard.right_shift)

"The character typed by `key`, for the small set allowed in labels, or `nothing`."
function keychar(key::Keyboard.Button, upper::Bool)
    name = string(key)
    if length(name) == 1 && 'a' ≤ name[1] ≤ 'z'
        upper ? uppercase(name[1]) : name[1]
    elseif length(name) == 2 && name[1] == '_' && isdigit(name[2])
        name[2]
    elseif key == Keyboard.minus
        '-'
    elseif key == Keyboard.period
        '.'
    elseif key == Keyboard.space
        ' '
    else
        nothing
    end
end

"Close every window once the event that asked for it has been handled."
quit!() = @async (sleep(0.05); GLMakie.closeall())

"""
Draw the stack. Traces are a neutral colour by default so that the multiplet colours stand
out; pass a `colormap` to colour them by plane instead.
"""
function draw!(ax::Axis, s::SpectrumSeries, v::StackView; colormap=nothing,
               color=(:black, 0.6))
    n = length(s)
    x = shifts(s)
    for i in eachindex(s)
        style = isnothing(colormap) ? (; color) :
                (; color=i, colormap, colorrange=(1, max(n, 2)))
        trace = lift(v.offsets, v.gain, v.revision) do offsets, g, _
            return Point2f.(x .+ offsets[i],
                            ypos.(v.spacing, g, v.norm, i, intensities(s, i)))
        end
        lines!(ax, trace; linewidth=1, style...)
    end
    # limits given high-to-low: Makie takes the axis direction from their order,
    # overriding xreversed, so this is what makes chemical shift run right to left
    lo, hi = extrema(x)
    limits!(ax, hi, lo, -0.5v.spacing, (n + 1) * v.spacing)
    return ax
end

"""
    openwindow(s, processing; labels, reference, refshift, coordinates, peaklist, radius, snr,
           gain, colormap, size, call) -> PeakTracker

Open the picking window on `s`, wait until it is closed, and return the analysis.
"""
function openwindow(s::SpectrumSeries, processing; labels=String[], reference="reference",
                    refshift=0.0, coordinates=nothing, peaklist=nothing, radius=8, snr=5,
                    gain=2, colormap=nothing, size=(1200, 850), call=nothing,
                    outputdir="out")
    coords = coordinatetable(coordinates, length(s))
    t = initialtracks(s, peaklist, cleanlabel(reference),
                      String[cleanlabel(l) for l in labels])

    GLMakie.activate!()   # CairoMakie, loaded for saving figures, would otherwise take over
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xreversed=true, xlabel=shiftlabel(s.nucleus), title=s.title,
              yticksvisible=false, yticklabelsvisible=false,
              xgridvisible=false, ygridvisible=false,
              leftspinevisible=false, rightspinevisible=false, topspinevisible=false,
              yzoomlock=true, ypanlock=true)   # the mouse zooms and pans chemical shift only
    for (row, text) in enumerate(HELP)
        Label(fig[row + 1, 1:2], text; fontsize=12, color=:grey40, tellwidth=false)
    end
    v = StackView(s; gain)
    draw!(ax, s, v; colormap)

    p = PeakTracker(s, coords, Observable(t), processing, Float64(refshift),
                    Observable(false), v, fig, ax, call, Observable(String(outputdir)))
    controls = GridLayout(fig[4, 1:2]; tellwidth=false)
    isnothing(processing) || broadeningcontrol!(controls[1, 1], p)
    textbox = outputcontrols!(controls[1, 2], p)
    attach!(p, textbox, fig[1, 2]; radius, snr)
    colsize!(fig.layout, 2, Relative(0.28))

    display(fig)
    while isopen(fig.scene)
        sleep(0.1)
    end
    GLMakie.closeall()
    return p
end

"""
    summarypanel!(position, p, current)

The panel beside the stack: the shifts (above) and heights (below) of the current line's
multiplet against the first coordinate, following the selection and every edit. See
[`summaryseries`](@ref).
"""
function summarypanel!(position, p::PeakTracker, current)
    grid = GridLayout(position)
    data = lift((_, c) -> summaryseries(p, c), p.tracks, current)
    colour = @lift(multipletcolour(max($data.m, REF)))
    axshift = Axis(grid[1, 1]; title=@lift($data.label), ylabel=@lift($data.ylabel),
                   xticklabelsvisible=false)
    axheight = Axis(grid[2, 1]; xlabel=@lift($data.xlabel), ylabel="height")
    linkxaxes!(axshift, axheight)
    scatter!(axshift, @lift($data.lines); color=@lift(($colour, 0.5)), markersize=6)
    scatterlines!(axshift, @lift($data.centre); color=colour, markersize=8)
    scatterlines!(axheight, @lift($data.heights); color=colour, markersize=6)
    # autolimits on series with no finite point would collapse the axes
    finite(points) = any(q -> all(isfinite, q), points)
    on(data; update=true) do d
        finite(d.centre) || finite(d.lines) || return
        autolimits!(axshift)
        finite(d.heights) && autolimits!(axheight)
        return
    end
    return grid
end

"Multiplets to start from: named ones, a saved peak list, or an earlier analysis's."
function initialtracks(s, ::Nothing, reference, labels)
    return Tracks([reference; filter(!=(reference), unique(labels))])
end

function initialtracks(s, path::AbstractString, reference, labels)
    return addlabels!(restore(s, path; reference), labels)
end

function initialtracks(s, previous::PeakTracker, reference, labels)
    t = previous.series === s ? copy(tracks(previous)) :
        restore(s, peaklisttable(previous)...; reference, path="the earlier analysis")
    return addlabels!(t, labels)
end

"Add an empty multiplet for each of `labels` not already present."
function addlabels!(t::Tracks, labels)
    for l in labels
        l in multipletlabels(t) || addmultiplet!(t, l)
    end
    return t
end

"""
A slider for the exponential line broadening. Changing it reprocesses every experiment;
picks keep their shifts, and their heights are read from the new data.
"""
function broadeningcontrol!(position, p::PeakTracker)
    proc = p.processing
    lbs = sort!(unique([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0,
                        10.0,
                        15.0, 20.0, proc.lb[]]))
    grid = SliderGrid(position,
                      (label="exponential line broadening", range=lbs, startvalue=proc.lb[],
                       format=x -> "$(x) Hz"); tellwidth=false, width=500)
    on(grid.sliders[1].value) do lb
        lb == proc.lb[] && return
        proc.lb[] = lb
        process!(p.series, proc, lb)
        notify(p.view.revision)
        reheight!(tracks(p), p.series)
        return notify(p.tracks)
    end
    return grid
end

"The output folder box and the Save, Load and Quit buttons; returns the box."
function outputcontrols!(position, p::PeakTracker)
    row = GridLayout(position)
    Label(row[1, 1], "output folder")
    box = Textbox(row[1, 2]; width=120, stored_string=p.outputdir[])
    on(box.stored_string) do text
        return p.outputdir[] = isempty(strip(text)) ? "out" : String(strip(text))
    end
    # Makie's Textbox commits only on enter; clicking away would otherwise lose the text
    on(box.focused) do focused
        return focused || (box.stored_string[] = box.displayed_string[])
    end
    savebutton = Button(row[1, 3]; label="Save")
    on(_ -> saveanalysis(p, joinpath(pwd(), p.outputdir[])), savebutton.clicks)
    loadbutton = Button(row[1, 4]; label="Load")
    on(loadbutton.clicks) do _
        path = joinpath(pwd(), p.outputdir[], "peaklist.csv")
        if !isfile(path)
            @warn "There is no peak list to load" path
            return
        end
        p.tracks[] = restore(p.series, path; reference=multiplets(tracks(p))[REF].label)
        @info "Loaded $(npicks(tracks(p))) picks from $path"
    end
    quitbutton = Button(row[1, 5]; label="(Q)uit")
    on(_ -> quit!(), quitbutton.clicks)
    return box
end

"""
The line above the stack: every multiplet in its colour, the current one bracketed with the
current line's number where it has several, the vacant slot, and then whatever is being
typed or chosen, or the brush settings.
"""
function statusline(t::Tracks, current, typed, joining, radius, aligned, snapping)
    parts = Any[]
    cm = multipletof(t, current)
    for (m, mp) in enumerate(multiplets(t))
        ls = navorder(t, m)
        name = length(ls) > 1 ? "$(mp.label) ×$(length(ls))" : mp.label
        if m == cm
            k = findfirst(==(current), ls)
            name = length(ls) > 1 ? "[$(mp.label) $k/$(length(ls))]" : "[$name]"
        end
        push!(parts, rich(name; color=multipletcolour(m), font=m == cm ? :bold : :regular),
              "   ")
    end
    vacant = current == 0
    push!(parts,
          rich(vacant ? "[+ $(vacantname(multipletlabels(t)))]" : "+";
               color=multipletcolour(nmultiplets(t) + 1),
               font=vacant ? :bold : :regular))
    if !isnothing(typed)
        push!(parts, "      rename to: $(typed)▏  (enter to confirm, esc to cancel)")
    elseif joining
        push!(parts,
              "      group: click a line to join it to $(t.multiplets[cm].label) " *
              "(esc to cancel)")
    else
        push!(parts,
              "      search ±$(round(Int, radius)) px   snap $(snapping ? "on" : "off")")
        aligned && push!(parts, "   aligned to $(t.multiplets[REF].label)")
    end
    return rich(parts...)
end

"""
    attach!(p, textbox, panel; radius=8, snr=5)

Wire up the mouse and keyboard, and draw the summary panel at `panel`. The current line is
the one that picking, moving, erasing, renaming and deleting act on; line `0` is the vacant
slot at the end, and picking there creates a new multiplet.
"""
function attach!(p::PeakTracker, textbox, panel; radius=8, snr=5)
    ax, v, s, found, aligned = p.axis, p.view, p.series, p.tracks, p.aligned
    n = length(s)

    current = Observable(isempty(linesof(found[], REF)) ? 0 : first(linesof(found[], REF)))
    renaming = Observable{Union{Nothing,String}}(nothing)
    grouping = Observable(false)
    target = Ref(p.refshift)
    r = Observable(Float32(radius))
    live = Observable(Point2f[])
    erasing = Observable(false)
    snapping = Observable(true)
    moving = Ref{Union{Nothing,Tuple{Int,Int}}}(nothing)
    history = Tracks[]
    # undo must not bring back heights read from the spectra as processed before
    on(_ -> foreach(h -> reheight!(h, s), history), v.revision)

    # every line in the order left and right step through them, then the vacant slot
    slots() = [[l for m in eachindex(multiplets(found[])) for l in navorder(found[], m)]; 0]
    function select!(step)
        list = slots()
        k = something(findfirst(==(current[]), list), length(list))
        return current[] = list[clamp(k + step, 1, length(list))]
    end
    # the current line, creating a multiplet for it first if the vacant slot is selected
    function active()
        if current[] == 0
            t = found[]
            current[] = addmultiplet!(t, vacantname(multipletlabels(t)))
        end
        return current[]
    end
    # after an undo or a load, the current line may no longer exist
    function validate!()
        return current[] == 0 || hasline(found[], current[]) || (current[] = first(slots()))
    end
    on(_ -> validate!(), found)
    checkpoint() = push!(history, copy(found[]))
    halfwidth(axis) = r[] / pixelscale(axis)[1]   # search radius in ppm
    # place exactly at the cursor? shift inverts the snapping setting
    precise(scene) = snapping[] == shifted(scene)
    place(scene, i, x, w) = precise(scene) ? exact(s, i, x) : snap(s, i, x, w)
    # the reference is frozen while the spectra are aligned to it
    frozen(l) = aligned[] && isreference(found[], l)
    function locked(l)
        frozen(l) || return false
        @warn "The reference can't be changed while the spectra are aligned to it; " *
              "press a to unalign first"
        return true
    end

    # keep the offsets in step with the reference picks
    function realign!()
        if aligned[]
            δ = referenceshifts(found[], n)
            if isnothing(δ)
                aligned[] = false   # re-enters realign! and clears the offsets
                return
            end
            offsets = target[] .- δ
            offsets == v.offsets[] || (v.offsets[] = offsets)   # redraw only if they move
        elseif any(!iszero, v.offsets[])
            v.offsets[] = zeros(n)
        end
    end
    on(_ -> realign!(), found)
    on(_ -> realign!(), aligned)
    function togglealignment!()
        if aligned[]
            aligned[] = false
            return
        end
        if isnothing(referenceshifts(found[], n))
            missingplanes = [i for i in 1:n if isnothing(centre(found[], REF, i))]
            @warn "Pick the reference in every plane before aligning" missingplanes
            return
        end
        return aligned[] = true
    end

    # overlays are created up front so that editing never triggers a limit reset
    layers = lift((t, c, g, o) -> overlay(t, c, g, o, v.spacing, v.norm),
                  found, current, v.gain, v.offsets)
    others = @lift($layers[1])
    siblings = @lift($layers[2])
    chosen = @lift($layers[3])
    bars = @lift($layers[4])
    summarypanel!(panel, p, current)
    highlight = lift(found, current) do t, c
        return multipletcolour(something(multipletof(t, c), nmultiplets(t) + 1))
    end
    linesegments!(ax, @lift($bars.segments); color=@lift(multipletcolour.($bars.colours)),
                  linewidth=1.5)
    linesegments!(ax, @lift($bars.path);
                  color=@lift(multipletcolour.($bars.pathcolours, 0.7)),
                  linewidth=1, linestyle=:dash)
    linesegments!(ax, @lift($others.segments);
                  color=@lift(multipletcolour.($others.segcolours, 0.6)), linewidth=1.5)
    scatter!(ax, @lift($others.points); color=@lift(multipletcolour.($others.colours, 0.6)),
             marker=:dtriangle, markersize=9)
    linesegments!(ax, @lift($siblings.segments); color=highlight, linewidth=2)
    scatter!(ax, @lift($siblings.points); color=highlight, marker=:dtriangle, markersize=11)
    linesegments!(ax, @lift($chosen.segments); color=highlight, linewidth=3)
    scatter!(ax, @lift($chosen.points); color=highlight, marker=:dtriangle, markersize=14,
             strokecolor=:black, strokewidth=1)
    lines!(ax, live; linewidth=@lift(2 * $r),
           color=@lift(($erasing ? Makie.RGBAf(0.5, 0.5, 0.5, 1) : $highlight, 0.35)))

    status = lift(statusline, found, current, renaming, grouping, r, aligned, snapping)
    ax.subtitle = status[]
    on(text -> ax.subtitle = text, status)

    # rename the current multiplet, or create one with this name at the vacant slot
    function rename!(name)
        t = found[]
        if current[] == 0
            name = cleanlabel(name)
            isempty(name) && return
            if name in multipletlabels(t)
                @warn "A multiplet called $name already exists"
                return
            end
            checkpoint()
            current[] = addmultiplet!(t, name)
        else
            m = multipletof(t, current[])
            checkpoint()
            if !relabel!(t, m, name)
                pop!(history)
                @warn "Can't rename to \"$name\": the name is blank or already taken"
            end
        end
        return notify(found)
    end

    function joinline!(k)
        grouping[] = false
        isnothing(k) && return
        t = found[]
        (locked(current[]) || locked(first(k))) && return
        multipletof(t, first(k)) == multipletof(t, current[]) && return
        checkpoint()
        group!(t, current[], first(k))
        return notify(found)
    end

    deactivate_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :picker) do event::MouseEvent, axis
        t = event.type
        scale = pixelscale(axis)
        if grouping[]
            t === MouseEventTypes.leftclick &&
                joinline!(hit(found[], v, event.data, scale, 10))
            return
        end
        if t === MouseEventTypes.leftdown
            moving[] = hit(found[], v, event.data, scale, 10, current[])
        elseif t === MouseEventTypes.leftdragstart
            if isnothing(moving[])
                erasing[] = shifted(axis.scene)
                live[] = [event.prev_data, event.data]
            elseif !locked(first(moving[]))
                checkpoint()
            end
        elseif t === MouseEventTypes.leftdrag
            if isnothing(moving[])
                push!(live[], event.data)
                notify(live)
            elseif !frozen(first(moving[]))
                (l, i) = moving[]
                setpick!(found[], l, exact(s, i, measured(v, i, event.data[1])))
                notify(found)
            end
        elseif t === MouseEventTypes.leftdragstop
            if isnothing(moving[])
                pts = live[]
                live[] = Point2f[]   # end the live stroke first, whatever happens next
                if erasing[]
                    if current[] != 0 && !locked(current[])
                        checkpoint()
                        band = densify(pts, r[] / 2, scale)
                        erase!(q -> covered(band, location(v, q), scale, r[]), found[],
                               current[])
                    end
                elseif !locked(current[])
                    picked = track(s, v, pts, r[], scale; snr)
                    if !isempty(picked)
                        checkpoint()
                        l = active()
                        foreach(q -> setpick!(found[], l, q), picked)
                    end
                end
                erasing[] = false
            else
                (l, i) = moving[]
                q = pickat(found[], l, i)
                frozen(l) || isnothing(q) ||
                    setpick!(found[], l, place(axis.scene, i, q.shift, halfwidth(axis)))
                moving[] = nothing
            end
            notify(found)
        elseif t === MouseEventTypes.leftclick
            moving[] = nothing
            locked(current[]) && return
            i = nearestrow(s, v, event.data, halfwidth(axis))
            isnothing(i) && return
            q = place(axis.scene, i, measured(v, i, event.data[1]), halfwidth(axis))
            isnothing(q) && return
            checkpoint()
            setpick!(found[], active(), q)
            notify(found)
        elseif t === MouseEventTypes.rightclick
            k = hit(found[], v, event.data, scale, 10, current[])
            (isnothing(k) || locked(first(k))) && return
            checkpoint()
            delete!(found[].picks, k)
            notify(found)
        end
    end

    on(events(ax.scene).keyboardbutton) do ev
        ev.action in (Keyboard.press, Keyboard.repeat) || return
        textbox.focused[] && return   # typing an output folder name

        # while renaming, keys type the new name and nothing else
        if !isnothing(renaming[])
            if ev.key in (Keyboard.enter, Keyboard.kp_enter)
                rename!(renaming[])
                renaming[] = nothing
            elseif ev.key == Keyboard.escape
                renaming[] = nothing
            elseif ev.key == Keyboard.backspace
                renaming[] = String(chop(renaming[]))
            else
                c = keychar(ev.key, shifted(ax.scene))
                isnothing(c) || (renaming[] = renaming[] * c)
            end
            return
        end
        if grouping[]
            ev.key == Keyboard.escape && (grouping[] = false)
            return
        end

        if ev.key == Keyboard.up
            v.gain[] *= STEP
        elseif ev.key == Keyboard.down
            v.gain[] /= STEP
        elseif ev.key == Keyboard.left
            select!(-1)
        elseif ev.key == Keyboard.right
            select!(1)
        end
        ev.action == Keyboard.press || return
        t = found[]
        if Keyboard._0 ≤ ev.key ≤ Keyboard._9
            m = Int(ev.key) - Int(Keyboard._0) + 1   # key 0 is the reference
            if m ≤ nmultiplets(t)
                current[] = first(navorder(t, m))
            elseif m == nmultiplets(t) + 1
                current[] = 0
            end
        elseif ev.key == Keyboard.r
            renaming[] = ""
        elseif ev.key == Keyboard.g
            current[] == 0 || locked(current[]) || (grouping[] = true)
        elseif ev.key == Keyboard.u
            l = current[]
            if l != 0 && !locked(l) && length(linesof(t, multipletof(t, l))) > 1
                checkpoint()
                ungroup!(t, current[])
                notify(found)
            end
        elseif ev.key == Keyboard.d
            if current[] != 0 && !locked(current[])
                checkpoint()
                deleteline!(t, current[])
                notify(found)
            end
        elseif ev.key == Keyboard.a
            togglealignment!()
        elseif ev.key == Keyboard.m
            snapping[] = !snapping[]
        elseif ev.key == Keyboard.left_bracket
            r[] = max(2.0f0, r[] - 2)
        elseif ev.key == Keyboard.right_bracket
            r[] += 2
        elseif ev.key == Keyboard.backspace && !isempty(history)
            found[] = pop!(history)
        elseif ev.key == Keyboard.s
            dir = mkpath(joinpath(pwd(), p.outputdir[]))
            path = joinpath(dir, "stack-$(Dates.format(now(), "yyyymmdd-HHMMSS")).pdf")
            save(path, stackplot(p); backend=CairoMakie, pt_per_unit=1)
            @info "Saved a figure of the displayed region" path = abspath(path)
        elseif ev.key == Keyboard.q
            quit!()
        end
    end

    return p
end
