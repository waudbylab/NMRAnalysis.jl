function process_mousebutton(expt, state, event)
    # consume input if in process of renaming a peak

    if state[:mode][] == :normal && event.button == Mouse.left &&
       event.action == Mouse.press
        pos = events(state[:gui][][:axcontour]).mouseposition[]
        pick_plt, pick_i = pick(state[:gui][][:axcontour], pos, 20)
        if pick_plt == state[:gui][][:pltinitialpeaks]
            @debug "event: initiate drag" pick_i
            state[:mode][] = :moving
            return Consume(true)
        end
    elseif state[:mode][] == :moving && event.button == Mouse.left &&
           event.action == Mouse.release
        @debug "event: stop dragging"
        state[:mode][] = :normal
        newpos = mouseposition(state[:gui][][:axcontour])
        movepeak!(expt, state[:current_peak_idx][], newpos)
        return Consume(true)
    end
    return Consume(false)
end

function process_mouseposition(expt, state, mousepos)
    if state[:mode][] == :moving
        newpos = mouseposition(state[:gui][][:axcontour])
        state[:initialpositions][][state[:current_peak_idx][]] = Point2f(newpos)
        notify(state[:initialpositions])
        return Consume()
    end
    return Consume(false)
end
