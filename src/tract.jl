"""
    tract()
    tract(trosy_filename, antitrosy_filename)

Analyze 1D TRACT relaxation experiments.

Interactively select integration and noise regions, then fit relaxation data to extract
relaxation rates and an effective correlation time.
    
- `trosy_filename`: Path to Bruker experiment folder for TROSY data
- `antitrosy_filename`: Path to Bruker experiment folder for anti-TROSY data

# Example

```julia
tract()
tract("path/to/trosy_experiment", "path/to/antitrosy_experiment")
```
"""
function tract()
    println("Current directory: $(pwd())")
    println()

    print("Enter path to TROSY experiment (i.e. Bruker experiment folder): ")
    trosy = readline()
    ispath(trosy) || throw(SystemError("No such file or directory"))

    print("Enter path to anti-TROSY experiment: ")
    ispath(trosy) || throw(SystemError("No such file or directory"))
    antitrosy = readline()

    return tract(trosy, antitrosy)
end

function tract(trosy, antitrosy)
    return tract(loadnmr(string(trosy)),
                 loadnmr(string(antitrosy)))
end

function tract(trosy::NMRData{T,2}, antitrosy::NMRData{T,2}) where {T}
    # numerical constants
    μ0 = 4π * 1e-7
    γH = gyromagneticratio(H1)
    γN = gyromagneticratio(N15)
    ħ = 6.626e-34 / 2π
    rNH = 1.02e-10
    ΔδN = 160e-6
    θ = 17 * π / 180

    # 1. get experiment parameters
    trosy = deepcopy(trosy) # work on copies of the data
    antitrosy = deepcopy(antitrosy)
    label!(trosy, "TROSY")
    label!(antitrosy, "Anti-TROSY")
    trosy = setrelaxtimes(trosy, acqus(trosy, :vdlist), "s")
    antitrosy = setrelaxtimes(antitrosy, acqus(antitrosy, :vdlist), "s")

    trosytau = data(trosy, 2)
    antitrosytau = data(antitrosy, 2)

    # calculate field-dependent quantities
    B0 = 2π * acqus(trosy, :bf1) / γH
    ωN = 2π * acqus(trosy, :bf3)

    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = B0 * γN * ΔδN / (3 * sqrt(2))
    f = p * c * (3cos(θ)^2 - 1)

    # 2. pick integration region

    plt = plot([trosy[:, 1],
                antitrosy[:, 1]]; grid=true)
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
    p2 = plot(trosy[:, 1]; linecolor=:black)
    hline!(p2, [0]; c=:grey, primary=false)
    plot!(p2, trosy[roi, 1]; fill=(0, :dodgerblue), linecolor=:navy,
          label="integration region")
    plot!(p2, trosy[noiseroi, 1]; fill=(0, :orange), linecolor=:red, label="noise region",
          legend=:topright)
    title!(p2, "Integration regions")
    display(p2)
    print("Displaying integration and noise regions. Press enter to continue.")
    readline()

    # 4. integrate regions
    trosynoise = vec(data(sum(trosy[noiseroi, :]; dims=F1Dim)))
    antitrosynoise = vec(data(sum(antitrosy[noiseroi, :]; dims=F1Dim)))
    trosynoise = std(trosynoise)
    antitrosynoise = std(antitrosynoise)

    trosyintegrals = vec(data(sum(trosy[roi, :]; dims=F1Dim)))
    antitrosyintegrals = vec(data(sum(antitrosy[roi, :]; dims=F1Dim)))

    # normalise by max value
    trosynoise /= maximum(trosyintegrals)
    antitrosynoise /= maximum(antitrosyintegrals)

    trosyintegrals /= maximum(trosyintegrals)
    antitrosyintegrals /= maximum(antitrosyintegrals)

    # 5. fit region
    model(t, p) = p[1] * exp.(-p[2] * t)
    p0 = [1.0, 5.0] # initial guesses for amplitude, R2
    trosyfit = curve_fit(model, trosytau, trosyintegrals, p0)
    antitrosyfit = curve_fit(model, antitrosytau, antitrosyintegrals, p0)
    trosypars = coef(trosyfit) .± stderror(trosyfit)
    antitrosypars = coef(antitrosyfit) .± stderror(antitrosyfit)
    trosyR2 = trosypars[2]
    antitrosyR2 = antitrosypars[2]
    trosyI0 = coef(trosyfit)[1]
    antitrosyI0 = coef(antitrosyfit)[1]

    # 6. calculate τc
    ΔR = antitrosyR2 - trosyR2
    ηxy = ΔR / 2

    # based on analytical solution in Mathematica
    # Solve[f (4*4/10*tc + (3*4/10*tc/(1 + (tc*\[Omega]N)^2))) == \[Eta], tc]
    x2 = 21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy^2 * ωN^8 + 625 * f^2 * ηxy^4 * ωN^10
    x = sqrt(x2)
    y3 = 1800 * f^2 * ηxy * ωN^4 + 125 * ηxy^3 * ωN^6 + 24 * sqrt(3) * x
    y = cbrt(y3)
    τc = (5 * ηxy) / (24 * f) -
         (336 * f^2 * ωN^2 - 25 * ηxy^2 * ωN^4) / (24 * f * ωN^2 * y) + y / (24 * f * ωN^2)
    τc = 1e9 * τc

    # 7. plot fitted region
    maxt = max(maximum(trosytau), maximum(antitrosytau))
    t = LinRange(0, maxt * 1.1, 100)
    yfit1 = model(t, coef(trosyfit))
    yfit2 = model(t, coef(antitrosyfit))

    p1 = scatter(1000trosytau, (trosyintegrals .± trosynoise) / trosyI0; label="TROSY",
                 frame=:box,
                 xlabel="Relaxation time / ms",
                 ylabel="Integrated signal",
                 title="TRACT: τc = $τc ns",
                 grid=nothing)
    scatter!(p1, 1000antitrosytau, (antitrosyintegrals .± antitrosynoise) / antitrosyI0;
             label="Anti-TROSY")
    plot!(p1, 1000t, yfit1 / trosyI0;
          label="TROSY fit (R₂ = $trosyR2 s⁻¹)",
          c=1,
          z_order=:back)
    plot!(p1, 1000t, yfit2 / antitrosyI0;
          label="Anti-TROSY fit (R₂ = $antitrosyR2 s⁻¹)",
          c=2,
          z_order=:back)

    println()
    @info """TRACT results

Current directory: $(pwd())
TROSY experiment: $(trosy[:filename])
Anti-TROSY experiment: $(antitrosy[:filename])

Integration region: $ppm1 - $ppm2 ppm
Noise region: $noise1 - $noise2 ppm

Fitted TROSY relaxation rate: $trosyR2 ns
Fitted anti-TROSY relaxation rate: $antitrosyR2 ns

Estimated τc: $τc ns

N.B. TRACT analysis assumes the protein is perfectly rigid. In the presence of flexibility or disorder, reported τc values will be underestimates.
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