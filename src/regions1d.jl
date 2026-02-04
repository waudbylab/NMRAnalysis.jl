"""
    get1dregionandnoise(spec::NMRData{T,1}) where {T}

Interactively select an integration region and a noise region from a 1D NMR spectrum (or slice).
The function will display the spectrum, prompt the user to enter the chemical shift values for 
the integration region and the noise region, and then return the corresponding selectors for
those regions.
"""
function get1dregionandnoise(spec::NMRData{T,1}) where {T}
    # show spectrum
    plt = plot(spec; grid=true)
    hline!(plt, [0]; c=:grey)
    display(plt)

    # prompt for integration region
    print("Defining integration region - please enter first chemical shift: ")
    ppm1 = readline()
    ppm1 = tryparse(Float64, ppm1)
    print("Defining integration region - please enter second chemical shift: ")
    ppm2 = readline()
    ppm2 = tryparse(Float64, ppm2)
    if ppm1 === nothing || ppm2 === nothing
        throw(ArgumentError("Invalid chemical shift value(s) entered for integration region. Please enter valid numbers."))
    end
    ppm1, ppm2 = minmax(ppm1, ppm2)

    # display integration region
    vspan!(plt, [ppm1, ppm2]; label="integration region", alpha=0.2)
    display(plt)

    # prompt for noise region
    print("Enter a chemical shift in the center of the noise region: ")
    noiseppm = readline()
    noiseppm = tryparse(Float64, noiseppm)
    if noiseppm === nothing
        throw(ArgumentError("Invalid chemical shift for noise region"))
    end

    # create integration region and noise selectors
    ppmrange = ppm2 - ppm1
    noise1 = noiseppm - 0.5 * ppmrange
    noise2 = noiseppm + 0.5 * ppmrange

    roi = ppm1 .. ppm2
    noiseroi = noise1 .. noise2

    # plot regions
    p2 = plot(spec[:, 1]; linecolor=:black)
    hline!(p2, [0]; c=:grey, primary=false)
    plot!(p2, spec[roi, 1]; fill=(0, :dodgerblue), linecolor=:navy,
          label="integration region")
    plot!(p2, spec[noiseroi, 1]; fill=(0, :orange), linecolor=:red, label="noise region",
          legend=:topright)
    title!(p2, "Integration regions")
    display(p2)
    print("Displaying integration and noise regions. Press enter to continue.")
    readline()

    # return selectors for regions
    return roi, noiseroi
end