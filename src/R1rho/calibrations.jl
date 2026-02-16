function setupR1rhopowers(calibration_experiment_file="")
    # ANSI escape code for magenta
    magenta = "\033[35m"
    reset = "\033[0m"

    # Get the terminal width
    term_width = displaysize(stdout)[2]
    line_break = repeat("-", term_width)

    if calibration_experiment_file == ""
        # Prompt the user for inputs with magenta color and caret
        println()
        println("Enter p1 (19F hard pulse power, in us):")
        print("> ")
        p1 = parse(Float64, readline()) * 1e-6 # convert to seconds

        println()
        println("Enter pldB1 (19F hard pulse power, in dB):")
        print("> ")
        pldb1 = parse(Float64, readline())
        pl1 = Power(pldb1, :dB)
    else
        calibration = analyse(calibration_experiment_file)
        p1 = Measurements.value(calibration.pulse90)
        pl1 = calibration.power_level
    end

    println()
    println("Input a list of spinlock strengths (in Hz) separated by commas, or press ENTER for a default list (100-15000 Hz):")
    print("> ")
    input = readline()

    if input == ""
        target_spinlock_strengths = [300, 500, 750, 1000, 1500, 2000, 3000, 4000,
                                     5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
                                     13000, 14000, 15000]
        println("Using default spinlock strengths (in Hz):\n[300, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]")
    else
        target_spinlock_strengths = parse.(Float64, strip.(split(input, ",")))
    end

    # Check for high spinlock strengths (above 10 kHz)
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

    # Calculate the final powers
    final_powers = convert_Hz_to_dB.(target_spinlock_strengths, pl1, p1)

    # Shuffle the final powers list and the corresponding spinlock strengths
    shuffled_indices = shuffle(1:length(final_powers))
    shuffled_final_powers = final_powers[shuffled_indices]
    shuffled_spinlock_strengths = target_spinlock_strengths[shuffled_indices]

    # Print the final powers in the specified format
    println()
    println("The list corresponds to the following spinlock strengths (Hz):\n",
            shuffled_spinlock_strengths)
    println()
    println("Copy & paste the list provided between the dashed lines.")
    println(line_break)
    println("dB")
    for power in shuffled_final_powers
        println(@sprintf("%.2f", power))
    end
    return println(line_break)
end