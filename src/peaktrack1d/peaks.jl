# Finding peaks in one plane: snapping to the tallest point nearby, and the position of a
# maximum between points.

"Indices of the points within `w` ppm of `x`."
window(s::SpectrumSeries, x, w) = findall(d -> abs(d - x) ≤ w, shifts(s))

"Walk uphill from index `k` to the nearest local maximum of `y`."
function climb(y, k)
    while k > firstindex(y) && y[k - 1] > y[k]
        k -= 1
    end
    while k < lastindex(y) && y[k + 1] > y[k]
        k += 1
    end
    return k
end

"""
    vertex(x, y, k) -> (shift, height)

The maximum near index `k` by parabolic interpolation through it and its two neighbours.
"""
function vertex(x, y, k)
    (k == firstindex(y) || k == lastindex(y)) && return x[k], y[k]
    a, b, c = y[k - 1], y[k], y[k + 1]
    d = (a - c) / (2 * (a - 2b + c))
    return x[k] + d * (x[k + 1] - x[k - 1]) / 2, b - (a - c) * d / 4
end

"Plane `i` picked at the tallest peak within `w` ppm of `x`, or `nothing` if the window is empty."
function snap(s::SpectrumSeries, i, x, w)
    idx = window(s, x, w)
    isempty(idx) && return nothing
    y = intensities(s, i)
    δ, h = vertex(shifts(s), y, climb(y, idx[argmax(y[idx])]))
    return Pick(i, δ, h)
end

"Plane `i` picked exactly at shift `x`, with the height of the nearest point."
function exact(s::SpectrumSeries, i, x)
    k = argmin(abs.(shifts(s) .- x))
    return Pick(i, x, intensities(s, i)[k])
end

"Signal-to-noise ratio of a pick."
snr(s::SpectrumSeries, q::Pick) = q.height / noiselevel(s, q.plane)
