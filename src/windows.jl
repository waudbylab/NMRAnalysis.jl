# Shared conventions for the interactive analysis windows.

"""
    windowname(analysis) -> String

The title of an analysis window, naming the package and the analysis:
`windowname("1D peak tracking") == "NMRAnalysis.jl - 1D peak tracking"`.
"""
windowname(analysis) = "NMRAnalysis.jl - $analysis"

"A name as it reads mid-sentence: a lower-case first letter, unless it begins an acronym."
function midsentence(name)
    name = String(name)
    return length(name) > 1 && isuppercase(name[2]) ? name : lowercasefirst(name)
end

"""
Where [`screenshot`](@ref) is saving the next analysis window, and how to show it, as
`(; file, peak, zoom, contour, xlims)`, or `nothing` while windows open as usual.
"""
const SCREENSHOT = Ref{Any}(nothing)

"""
    screenshot(f, file; peak=nothing, zoom=1, contour=0, xlims=nothing) -> f()

Run `f`, saving the analysis window it opens to the image `file` instead of opening it. The
window is drawn once the fit of the peaks or regions it starts with has finished, and the
analysis then returns as if the window had been closed straight away. The documentation's
screenshots are made this way.

A 2D window shows the peak labelled `peak`, or the first peak, magnified `zoom` times about
it, with the contours raised (positive) or lowered (negative) by `contour` steps of the ↑
and ↓ keys. A 1D window shows the chemical shift range `xlims`, in ppm, as `(lo, hi)`.

# Example
```julia
NMRAnalysis.screenshot("relaxation1d.png"; xlims=(0, 5)) do
    reopen(exampleanalysis("relaxation1d_aas"))
end
```
"""
function screenshot(f, file::AbstractString; peak=nothing, zoom=1, contour=0,
                    xlims=nothing)
    previous = SCREENSHOT[]
    SCREENSHOT[] = (; file=abspath(file), peak, zoom, contour, xlims)
    try
        return f()
    finally
        SCREENSHOT[] = previous
    end
end

"""
    capturewindow(fig; ready=() -> nothing) -> Bool

Where [`screenshot`](@ref) asked for it, run `ready()`, save `fig` and return true, so the
caller returns without opening the window. Otherwise false.
"""
function capturewindow(fig; ready=() -> nothing)
    shot = SCREENSHOT[]
    isnothing(shot) && return false
    ready()
    mkpath(dirname(shot.file))
    GLMakie.save(shot.file, fig; backend=GLMakie, px_per_unit=1)
    return true
end
