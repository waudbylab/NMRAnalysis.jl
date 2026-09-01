"""
    diffusion1d(coherence=SQ(H1))
    diffusion1d(filename, coherence=SQ(H1))
    diffusion1d(nmrdata, coherence=SQ(H1))

Analyze NMR diffusion experiments by fitting signal decay to extract diffusion coefficients.

Interactively configure gradient parameters, select integration and noise regions, then fit
diffusion data to the Stejskal-Tanner equation to extract diffusion coefficients and hydrodynamic radii.

- `filename`: Path to Bruker experiment folder
- `nmrdata`: NMRData object (2D: chemical shift × gradient strength)
- `coherence`: Coherence type for gyromagnetic ratio calculation (default: `SQ(H1)`)

# Example

```julia
diffusion1d()
diffusion1d("path/to/experiment")
```
"""
function diffusion1d(coherence=SQ(H1))
    println("Current directory: $(pwd())")
    println()

    print("Enter path to diffusion experiment (i.e. Bruker experiment folder): ")
    experiment = readline()
    ispath(experiment) || throw(SystemError("No such file or directory"))

    return diffusion1d(experiment, coherence)
end

function diffusion1d(filename, coherence=SQ(H1))
    return diffusion1d(loadnmr(string(filename)), coherence)
end

function diffusion1d(spec::NMRData{T,2}, coherence=SQ(H1)) where {T}
    spec = deepcopy(spec) # work on copy of the data
    label!(spec, "Diffusion")

    # 1. get experiment parameters
    td = spec[2, :npoints]

    println("Parsing experiment parameters...")
    γ = gyromagneticratio(coherence)
    δ = acqus(spec, :p, 30) * 2.0   # gradient pulse length
    Δ = acqus(spec, :d, 20)         # diffusion delay
    gpnam = acqus(spec, :gpnam, 6)  # try to identify shape factor
    σ = 1
    if length(gpnam) ≥ 4
        if gpnam[1:4] == "SMSQ"
            σ = 0.9
        elseif gpnam[1:4] == "SINE"
            σ = 0.6366
        end
    end
    temp = acqus(spec, :te)
    solvent = acqus(spec, :solvent)
    if solvent == "D2O"
        solvent = :d2o
    elseif solvent == "H2O+D2O"
        solvent = :h2o
    else
        solvent = missing
    end

    print("Gradient pulse length δ = $(1e6*δ) μs (2*p30). Press enter to confirm or type correct value (in μs): ")
    response = readline()
    if length(response) > 0
        δ = tryparse(Float64, response) * 1e-6
    end

    print("Diffusion delay Δ = $Δ s (d20). Press enter to confirm or type correct value (in s): ")
    response = readline()
    if length(response) > 0
        Δ = tryparse(Float64, response)
    end

    print("Gradient shape factor σ = $σ (gpnam6 = $gpnam). Press enter to confirm or type correct value: ")
    response = readline()
    if length(response) > 0
        σ = tryparse(Float64, response)
    end

    Gmax = 0.55
    print("Max. gradient strength Gmax = $Gmax T m⁻¹ (typical value for Bruker systems). Press enter to confirm or type correct value (in T m⁻¹): ")
    response = readline()
    if length(response) > 0
        Gmax = tryparse(Float64, response)
    end

    print("Enter initial gradient strength (%): ")
    response = readline()
    g1 = tryparse(Float64, response) * 0.01

    print("Enter final gradient strength (%): ")
    response = readline()
    g2 = tryparse(Float64, response) * 0.01

    print("Enter gradient ramp type ('l' / 'q' / 'e'): ")
    response = readline()
    length(response) == 1 || throw(ArgumentError("Invalid input"))
    response = lowercase(response)[1]
    if response == 'l'
        g = LinRange(g1, g2, td)
    elseif response == 'q'
        throw(ArgumentError("Unsupported ramp shape"))
    elseif response == 'e'
        throw(ArgumentError("Unsupported ramp shape"))
    else
        throw(ArgumentError("Invalid input"))
    end

    spec = setgradientlist(spec, g, Gmax)

    # 2. pick integration region
    plt = plot(spec[:, 1]; grid=true)
    hline!(plt, [0]; c=:grey)
    display(plt)

    print("Defining integration region - please enter first chemical shift: ")
    ppm1 = readline()
    ppm1 = tryparse(Float64, ppm1)
    print("Defining integration region - please enter second chemical shift: ")
    ppm2 = readline()
    ppm2 = tryparse(Float64, ppm2)

    # 3. pick noise region
    vspan!(plt, [ppm1, ppm2]; label="integration region", alpha=0.2)
    display(plt)

    print("Enter a chemical shift in the center of the noise region: ")
    noiseppm = readline()
    noiseppm = tryparse(Float64, noiseppm)

    # create integration region and noise selectors
    ppm1, ppm2 = minmax(ppm1, ppm2)

    ppmrange = ppm2 - ppm1
    noise1 = noiseppm - 0.5ppmrange
    noise2 = noiseppm + 0.5ppmrange

    roi = ppm1 .. ppm2
    noiseroi = noise1 .. noise2

    # plot region
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

    # 4. integrate regions
    noise = vec(data(sum(spec[noiseroi, :]; dims=F1Dim)))
    noise = std(noise)

    integrals = vec(data(sum(spec[roi, :]; dims=F1Dim)))

    # normalise by max value
    noise /= maximum(integrals)
    integrals /= maximum(integrals)

    # 5. fit region
    model(g, p) = p[1] * exp.(-(γ * δ * σ * g * Gmax) .^ 2 .* (Δ - δ / 3) .* p[2] .* 1e-10)

    p0 = [1.0, 1.0] # initial guesses for amplitude, D
    fit = curve_fit(model, g, integrals, p0)
    pars = coef(fit) .± stderror(fit)
    I0 = coef(fit)[1]
    D = pars[2] * 1e-10

    # 7. plot fitted region
    maxg = maximum(g)
    x = LinRange(0, maxg * 1.1, 100)
    yfit = model(x, coef(fit))

    p1 = scatter(g * Gmax, (integrals .± noise) / I0; label="observed",
                 frame=:box,
                 xlabel="Gradient strength / T m⁻¹",
                 ylabel="Integrated signal",
                 title="Diffusion: D = $D m² s⁻¹",
                 grid=nothing)
    plot!(p1, x * Gmax, yfit / I0;
          label="fit",
          z_order=:back)

    # 8. estimate rH (if solvent permits)
    if !ismissing(solvent)
        η = viscosity(solvent, temp)
        kB = 1.38e-23
        rH = kB * temp / (6π * η * 0.001 * D) * 1e10 # in Å
    else
        η = missing
        rH = missing
    end

    println()
    @info """diffusion results

Current directory: $(pwd())
Experiment: $(spec[:filename])

Integration region: $ppm1 - $ppm2 ppm
Noise region: $noise1 - $noise2 ppm

Solvent: $solvent
Temperature: $temp K
Expected viscosity: $η mPa s

Fitted diffusion coefficient: $D m² s⁻¹
Calculated hydrodynamic radius: $rH Å
"""
    display(p1)

    # 8. make final plot

    println()
    print("Enter a filename to save figure (press enter to skip): ")
    outputfilename = readline()
    if length(outputfilename) > 0
        plot!(p1; title="")
        savefig(p1, outputfilename)
        println("Figure saved to $outputfilename.")
    end

    return p1
end