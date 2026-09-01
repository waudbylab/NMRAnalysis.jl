"extend LsqFit.curve_fit to accept ydata with uncertainties"
function LsqFit.curve_fit(model,
                          xdata::AbstractArray,
                          ydata::AbstractArray{Measurement{T}} where {T},
                          p0::AbstractArray;
                          inplace=false,
                          kwargs...,)
    y = Measurements.value.(ydata)
    ye = Measurements.uncertainty.(ydata)
    wt = ye .^ -2
    return curve_fit(model, xdata, y, wt, p0; inplace=inplace, kwargs...)
end
