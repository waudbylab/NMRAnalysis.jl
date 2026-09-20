"""
    setupR1rhopowers(calibration="")

Interactively calculate spin-lock power levels (in Watts) for an R1ρ relaxation
dispersion experiment, ready to paste into the acquisition software's power list.

- `calibration`: one or more 1D nutation calibration experiments, given as paths or
  experiment numbers, or a ready-made [`B1Calibration`](@ref). Each experiment is fitted in
  the analysis window, so the fits can be checked before powers are computed from them;
  several of them, recorded at different power levels, give a calibration curve and the
  powers follow it rather than the ideal √W law. If omitted, you are prompted to enter p1
  and pldB1 manually and the ideal law is assumed.

You are then prompted for the target spin-lock strengths, either as an explicit
comma-separated list (in Hz) or, if left blank, a minimum/maximum bound used to filter a
built-in default list. Spin-lock strengths above 10 kHz trigger a confirmation prompt, to
guard against probe damage.

The resulting powers are shuffled (so that later fitting isn't biased by monotonically
ordered spin-lock strengths) and printed together with the corresponding shuffled list of
spin-lock strengths, ready to be copied into the experiment setup.

# Example
```julia
# read pulse calibration from an experiment, then prompt for spin-lock strengths
setupR1rhopowers("examples/calibration/1")

# calibrate the amplifier's response over several power levels first
setupR1rhopowers(["examples/calibration/1", "examples/calibration/2"])

# prompt for pulse parameters and spin-lock strengths interactively
setupR1rhopowers()
```
"""
function setupR1rhopowers(calibration="")
    # ANSI escape code for magenta
    magenta = "\033[35m"
    reset = "\033[0m"

    # Get the terminal width
    term_width = displaysize(stdout)[2]
    line_break = repeat("-", term_width)

    cal = b1calibration(calibration)
    println()
    println("Using calibration: ", cal)

    println()
    println("Input a list of spin-lock strengths (in Hz) separated by commas, or press ENTER for a default list:")
    print("> ")
    input = readline()

    if input == ""
        println("Input minimum spin-lock power (in Hz) [300 Hz]:")
        print("> ")
        input = readline()
        if input == ""
            min_spinlock_strength = 300
        else
            min_spinlock_strength = parse(Float64, input)
        end
        println("Input maximum spin-lock power (in Hz) [8000 Hz]:")
        print("> ")
        input = readline()
        if input == ""
            max_spinlock_strength = 8000
        else
            max_spinlock_strength = parse(Float64, input)
        end

        target_spinlock_strengths = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 4000,
                                     5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
                                     13000, 14000, 15000]
        target_spinlock_strengths = filter(x -> x >= min_spinlock_strength &&
                                                x <= max_spinlock_strength,
                                           target_spinlock_strengths)
        println("Using spin-lock strengths (in Hz):")
        println(target_spinlock_strengths)
    else
        target_spinlock_strengths = parse.(Float64, strip.(split(input, ",")))
    end

    # Check for high spin-lock strengths (above 10 kHz)
    max_spinlock_strength = maximum(target_spinlock_strengths)
    if max_spinlock_strength > 10000
        println()
        println("$(magenta)WARNING - high spin-lock powers may cause damage to your probe! 
Check the spin-lock duration is within acceptable power limits. 
Maximum spin-lock strength will be $max_spinlock_strength Hz. 
Type 'yes' to proceed. Do you want to proceed? (yes/no):$reset")
        println()
        print("> ")
        confirmation = readline()
        if lowercase(confirmation) != "yes"
            println("Operation cancelled.")
            return
        end
    end

    # Calculate the final powers, inverting the calibration curve
    final_powers_W = watts.(Power.(target_spinlock_strengths, cal))

    # Shuffle the final powers list and the corresponding spin-lock strengths
    shuffled_indices = shuffle(1:length(final_powers_W))
    shuffled_final_powers = final_powers_W[shuffled_indices]
    shuffled_spinlock_strengths = target_spinlock_strengths[shuffled_indices]

    # Print the final powers in the specified format
    println()
    println("The list corresponds to the following spin-lock strengths (Hz):\n",
            shuffled_spinlock_strengths)
    println()
    println("Copy & paste the list provided between the dashed lines.")
    println(line_break)
    println("Watt")
    for power in shuffled_final_powers
        println(@sprintf("%.10f", power))
    end
    return println(line_break)
end

"""
    b1calibration(x) -> B1Calibration

The calibration `setupR1rhopowers` works from: one ready-made, or the nutation calibration
experiment(s) to measure one from, or `""` to ask for the hard pulse and its power level
and assume the ideal √W power law.

Each calibration experiment is fitted in the analysis window rather than headlessly: the
powers about to be pasted into the spectrometer are only as good as these fits, so they are
worth looking at.
"""
b1calibration(cal::B1Calibration) = cal
function b1calibration(x::AbstractVector)
    return B1Calibration(calibration1d(x); source=join(string.(x), ", "))
end

function b1calibration(x)
    x == "" || return b1calibration([x])
    println()
    println("Enter p1 (hard pulse length, in us):")
    print("> ")
    p1 = parse(Float64, readline()) * 1e-6 # convert to seconds

    println()
    println("Enter pldB1 (hard pulse power, in dB):")
    print("> ")
    pldb1 = parse(Float64, readline())
    # A 90° pulse of length p1 is a field of 1/(4·p1); with one measurement the power law
    # is assumed ideal, which is what this routine has always done.
    return B1Calibration([Power(pldb1, :dB)], [1 / (4p1)]; source="p1/pldB1 as entered")
end
