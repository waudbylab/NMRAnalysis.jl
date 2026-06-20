function process_keyboardbutton(expt, state, event)
    @debug "keyboard event: $event"
    g = state[:gui][]
    if state[:mode][] == :normal && event.action == Keyboard.press
        if ispressed(g[:fig], Keyboard.a)
            pos = mouseposition(g[:axcontour])
            addpeak!(expt, Point2f(pos))
            state[:current_peak_idx][] = length(expt.peaks[]) # select the new peak
        elseif ispressed(g[:fig], Keyboard.d)
            idx = state[:current_peak_idx][]
            if idx > 0
                state[:current_peak_idx][] = 0
                deletepeak!(expt, idx)
            end
        elseif ispressed(g[:fig], Keyboard.r)
            if state[:current_peak_idx][] > 0
                renamepeak!(expt, state, :keyboard)
            end
        elseif ispressed(g[:fig], Keyboard.up)
            g[:basecontour][] *= g[:contourscale][]
        elseif ispressed(g[:fig], Keyboard.down)
            g[:basecontour][] /= g[:contourscale][]
        elseif ispressed(g[:fig], Keyboard.left)
            # left
            i = state[:current_slice][]
            if i > 1
                i -= 1
                set_close_to!(state[:gui][][:sliderslice], i)
            end
        elseif ispressed(g[:fig], Keyboard.right)
            # right
            i = state[:current_slice][]
            if i < nslices(expt)
                i += 1
                set_close_to!(state[:gui][][:sliderslice], i)
            end
        end
    elseif state[:mode][] == :renaming || state[:mode][] == :renamingstart
        if event.action == Keyboard.press && event.key == Keyboard.enter
            state[:current_peak][].label[] = state[:current_peak][].label[][1:(end - 1)]
            state[:mode][] = :normal
            notify(expt.peaks)
            return Consume()
        elseif event.action == Keyboard.press && event.key == Keyboard.backspace
            if length(state[:current_peak][].label[]) > 1
                state[:current_peak][].label[] = state[:current_peak][].label[][1:(end - 2)] * "‸"
                notify(expt.peaks)
                return Consume()
            end
        elseif event.action == Keyboard.press && event.key == Keyboard.escape
            # restore previous label
            state[:current_peak][].label[] = state[:oldlabel][]
            notify(expt.peaks)

            state[:mode][] = :normal
            return (Consume())
        end
        return Consume(false)
    elseif state[:mode][] == :fitting
        if event.action == Keyboard.press && event.key == Keyboard.escape
            # cancel the in-flight fit - the running fit's residual will see the bumped
            # generation and abort. No new fit will run, so clear the fitting status
            # here to restore the UI to :normal.
            state[:fit_generation][] += 1
            state[:mode][] = :normal
            return Consume()
        end
        return Consume(false)
    end
end

function process_unicode_input(expt, state, character)
    @debug "Processing unicode input: $character"
    if state[:mode][] == :renamingstart
        state[:mode][] = :renaming
        if character == 'r'
            # discard 'r' character hanging over from initial keypress
            return Consume()
        end
    end
    if state[:mode][] == :renaming
        state[:current_peak][].label[] = state[:current_peak][].label[][1:(end - 1)] *
                                         character * "‸"
        notify(expt.peaks)
        return Consume()
    end
    return Consume(false)
end