function gui!(state)
    state[:gui] = Dict{Symbol,Any}()
    gui = state[:gui]
    fig = Figure(; size=(1200, 900))

    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    c3 = Makie.wong_colors()[3]
    c4 = Makie.wong_colors()[4]
    c5 = Makie.wong_colors()[5]
    c6 = Makie.wong_colors()[6]
    c7 = Makie.wong_colors()[7]

    top_panel = fig[1, 1] = GridLayout()
    bottom_panel = fig[2, 1] = GridLayout()
    right_panel = fig[1:2, 2] = GridLayout()

    gui[:specplottitle] = lift(state[:currentseries]) do i
        return "Observed spectrum (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz) – drag to move peak/noise positions"
    end
    ax_spectra = Axis(top_panel[1, 1];
                      xreversed=true,
                      yzoomlock=true,
                      ypanlock=true,
                      xrectzoom=true,
                      yrectzoom=false,
                      xlabel="Chemical shift (ppm)",
                      ylabel="Intensity",
                      title=gui[:specplottitle])
    hideydecorations!(ax_spectra)

    gui[:noisespan] = lift(state[:noiseppm], state[:dx]) do noise, dx
        return noise - dx / 2, noise + dx / 2
    end
    gui[:peakspan] = lift(state[:peakppm], state[:dx]) do peak, dx
        return peak - dx / 2, peak + dx / 2
    end
    hlines!(ax_spectra, [0]; color=:grey)
    gui[:peakspan] = vspan!(ax_spectra, @lift($(gui[:peakspan])[1]),
                            @lift($(gui[:peakspan])[2]); alpha=0.7, label="Peak")
    gui[:noisespan] = vspan!(ax_spectra, @lift($(gui[:noisespan])[1]),
                             @lift($(gui[:noisespan])[2]); alpha=0.7, label="Noise",
                             color=c4)
    lines!(ax_spectra, state[:currentspectrum]; label="Observed")
    axislegend(ax_spectra; position=:lt)

    input_panel = right_panel[1, 1] = GridLayout()

    input_panel[1, 1] = Label(fig, "Series:")
    slider_current = input_panel[1, 2:3] = Slider(fig; range=1:state[:nseries], width=150)
    gui[:slider_current] = slider_current
    connect!(state[:currentseries], slider_current.value)

    input_panel[3, 1] = Label(fig, "Peak position (ppm):")
    text_peakppm = input_panel[3, 2:3] = Textbox(fig;
                                                 stored_string=string(round(state[:peakppm][];
                                                                            digits=2)),
                                                 validator=Float64, width=150,
                                                 textpadding=(4, 4, 4, 4))
    gui[:text_peakppm] = text_peakppm
    on(text_peakppm.stored_string) do s
        return state[:peakppm][] = parse(Float64, s)
    end

    input_panel[4, 1] = Label(fig, "Noise position (ppm):")
    text_noiseppm = input_panel[4, 2:3] = Textbox(fig;
                                                  stored_string=string(round(state[:noiseppm][];
                                                                             digits=2)),
                                                  validator=Float64, width=150,
                                                  textpadding=(4, 4, 4, 4))
    gui[:text_noiseppm] = text_noiseppm
    on(text_noiseppm.stored_string) do s
        return state[:noiseppm][] = parse(Float64, s)
    end

    input_panel[2, 1] = Label(fig, "Integration width (ppm):")
    text_dx = input_panel[2, 2] = Textbox(fig;
                                          stored_string=string(round(state[:dx][];
                                                                     digits=2)),
                                          validator=Float64, width=50, textpadding=(4, 4, 4, 4))
    gui[:text_dx] = text_dx
    on(text_dx.stored_string) do s
        return state[:dx][] = parse(Float64, s)
    end
    button_optimisewidth = input_panel[2, 3] = Button(fig; label="Optimise")
    on(button_optimisewidth.clicks) do _
        optimisewidth!(state)
        return gui[:text_dx].displayed_string[] = string(round(state[:dx][]; digits=3))
    end

    input_panel[5, 1] = Label(fig, "Initial I0:")
    text_I0 = input_panel[5, 2:3] = Textbox(fig;
                                            stored_string=string(round(state[:initialI0][];
                                                                       digits=1)),
                                            validator=Float64, width=150, textpadding=(4, 4, 4, 4))
    gui[:text_I0] = text_I0
    on(text_I0.stored_string) do s
        return state[:initialI0][] = parse(Float64, s)
    end
    on(state[:intensities]) do I
        state[:initialI0][] = I[1]
        return text_I0.displayed_string[] = string(round(I[1]; digits=1))
    end

    input_panel[6, 1] = Label(fig, "Initial R2,0 (s⁻¹):")
    text_R20 = input_panel[6, 2:3] = Textbox(fig;
                                             stored_string=string(round(state[:initialR20][];
                                                                        digits=1)),
                                             validator=Float64, width=150, textpadding=(4, 4, 4, 4))
    gui[:text_R20] = text_R20
    on(text_R20.stored_string) do s
        return state[:initialR20][] = parse(Float64, s)
    end

    input_panel[7, 1] = Label(fig, "Initial Rex (s⁻¹):")
    text_Rex = input_panel[7, 2:3] = Textbox(fig;
                                             stored_string=string(round(state[:initialRex][];
                                                                        digits=1)),
                                             validator=Float64, width=150, textpadding=(4, 4, 4, 4))
    gui[:text_Rex] = text_Rex
    on(text_Rex.stored_string) do s
        return state[:initialRex][] = parse(Float64, s)
    end

    input_panel[8, 1] = Label(fig, "Initial kex (s⁻¹):")
    text_kex = input_panel[8, 2:3] = Textbox(fig;
                                             stored_string=string(round(exp(state[:initiallnk][]);
                                                                        digits=1)),
                                             validator=Float64, width=150, textpadding=(4, 4, 4, 4))
    gui[:text_kex] = text_kex
    on(text_kex.stored_string) do s
        return state[:initiallnk][] = log(parse(Float64, s))
    end

    input_panel[9, 1] = Label(fig, "Δδ stdev (ppm):")
    text_σΔδ = input_panel[9, 2:3] = Textbox(fig;
                                             stored_string=string(round(state[:σΔδ][];
                                                                        digits=1)),
                                             validator=Float64, width=150, textpadding=(4, 4, 4, 4))
    gui[:text_σΔδ] = text_σΔδ
    on(text_σΔδ.stored_string) do s
        return state[:σΔδ][] = parse(Float64, s)
    end

    gui[:fitplottitle] = lift(state[:currentseries]) do i
        return "Peak integrals (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz)"
    end
    ax_fit = Axis(bottom_panel[1, 1];
                  xlabel="TSL (ms)",
                  ylabel="Peak integral",
                  title=gui[:fitplottitle])
    plt_obserr = errorbars!(ax_fit, state[:currenterror])
    plt_obsscat = scatter!(ax_fit, state[:currentscatter]; label="Observed")
    plt_glob = lines!(ax_fit, state[:currentfit]; label="Global fit", color=c1)
    plt_noex = lines!(ax_fit, state[:currentfit_null]; label="No exchange model", color=c3,
                      linestyle=:dash)
    plt_expfit = lines!(ax_fit, state[:currentexpfit]; label="Exponential fit", color=c2)
    plt_resierr = errorbars!(ax_fit, state[:residualerror]; color=c4)
    plt_resis = scatter!(ax_fit, state[:residualpoints]; color=c4, label="Residuals")
    plt_noexresierr = errorbars!(ax_fit, state[:residualerror_null]; color=c5)
    plt_noexresis = scatter!(ax_fit, state[:residualpoints_null]; color=c5,
                             label="Residuals (no ex)")

    # axislegend(ax_fit; position=:rt)
    Legend(bottom_panel[2, 1],
           ax_fit;
           nbanks=2, orientation=:horizontal)

    ax_fit_R1rho = Axis(bottom_panel[1, 2];
                        xlabel="νSL (kHz)",
                        ylabel="R1rho (s⁻¹)",
                        title="Dispersion curve")
    hlines!(ax_fit_R1rho, [0]; linewidth=0)
    errorbars!(ax_fit_R1rho, state[:expfiterror]; color=c2)
    scatter!(ax_fit_R1rho, state[:expfitpoints]; label="Exponential fits", color=c2)
    lines!(ax_fit_R1rho, state[:fitR1rho]; label="Global fit", color=c1)
    lines!(ax_fit_R1rho, state[:fitR1rho_null]; label="No exchange model", color=c3,
           linestyle=:dash)
    # axislegend(ax_fit_R1rho; position=:rt)
    Legend(bottom_panel[2, 2],
           ax_fit_R1rho;
           orientation=:horizontal)

    gui[:results_text] = lift(state[:ftest], state[:σΔδ]) do ftest, _
        I0 = state[:fitI0][]
        R20 = state[:fitR20][]
        Rex = state[:fitRex][]
        K = state[:fitK][]
        kex = state[:fitkex][]
        R20_null = state[:fitR20_null][]

        f_stat, p_value = ftest

        # Format full model results
        full_model = ["Exchange model results:",
                      "I₀ = $I0",
                      "R₂₀ = $R20 s⁻¹",
                      "Rₑₓ = $Rex s⁻¹",
                      "K = $K s⁻¹",
                      "kₑₓ = $kex s⁻¹"]

        # Format null model results
        null_model = ["No-exchange model results:",
                      "R₂₀ = $R20_null s⁻¹"]

        # Model comparison statistics
        comparison = ["Model comparison:",
                      "F-statistic = $(round(f_stat, digits=2)); p-value = $(round(p_value, digits=4))"]

        # Add significance interpretation
        if p_value < 0.05
            push!(comparison, "Exchange significant (p < 0.05)")
        else
            push!(comparison, "Exchange not significant (p ≥ 0.05)")
        end

        # Combine all sections with blank lines between them
        all_text = vcat(full_model, [""], null_model, [""], comparison)

        return join(all_text, "\n")
    end

    results_panel = right_panel[2, 1] = GridLayout(; tellheight=false)
    Label(results_panel[1, 1:2], gui[:results_text])
    # Label(results_panel[2,1:2], lift(x->"I0: $x", state[:fitI0]))
    # Label(results_panel[3,1:2], lift(x->"R2,0 (s⁻¹): $x", state[:fitR20]))
    # Label(results_panel[4,1:2], lift(x->"Rex (s⁻¹): $x", state[:fitRex]))
    # Label(results_panel[5,1:2], lift(x->"kex (s⁻¹): $(exp(x))", state[:fitlnk]))
    results_panel[2, 1:2] = Label(fig, "Working directory:\n$(pwd())")
    results_panel[3, 1] = Label(fig, "Output folder:")
    text_out = results_panel[3, 2] = Textbox(fig; stored_string="out", width=150,
                                             textpadding=(4, 4, 4, 4))
    gui[:text_out] = text_out
    on(text_out.stored_string) do s
        return state[:outputdir][] = s
    end
    button_save = results_panel[4, 1:2] = Button(fig; label="Save results")
    on(button_save.clicks) do _
        return savefig!(state)
    end

    gui[:dragging] = :nothing
    on(events(ax_spectra).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                if mouseover(fig, gui[:peakspan])
                    gui[:dragging] = :peak
                elseif mouseover(fig, gui[:noisespan])
                    gui[:dragging] = :noise
                else
                    gui[:dragging] = :nothing
                end
                return Consume(gui[:dragging] != :nothing)
            elseif event.action == Mouse.release
                # Exit dragging
                gui[:dragging] = :nothing
                return Consume(false)
            end
        end
    end
    on(events(fig).mouseposition; priority=2) do mp
        if gui[:dragging] != :nothing
            p = mouseposition(ax_spectra)
            if gui[:dragging] == :peak
                state[:peakppm][] = p[1]
                text_peakppm.displayed_string[] = string(round(p[1]; digits=3))
            elseif gui[:dragging] == :noise
                state[:noiseppm][] = p[1]
                text_noiseppm.displayed_string[] = string(round(p[1]; digits=3))
            end
            return Consume(true)
        end
        return Consume(false)
    end

    display(fig)
    optimisewidth!(state) # optimise width at start
    gui[:text_dx].displayed_string[] = string(round(state[:dx][]; digits=3))
    autolimits!(ax_fit)
    autolimits!(ax_fit_R1rho)

    while true #!state["should_close"][]
        sleep(0.1)
        if !isopen(fig.scene)
            break
        end
    end

    return GLMakie.closeall()
end

function savefig!(state)
    outputdir = joinpath(pwd(), state[:outputdir][])
    @info "Saving results to $outputdir"
    if !isdir(outputdir)
        mkdir(outputdir)
    else
        # move existing files to a backup folder
        backupdir = outputdir * "_previous"
        @info "Backing up previous results to $backupdir"
        if isdir(backupdir)
            rm(backupdir; recursive=true)
        end
        mv(outputdir, backupdir)
        mkdir(outputdir)
    end

    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]

    # dispersion fit
    fig = Figure()
    ax_fit_R1rho = Axis(fig[1, 1];
                        xlabel="νSL (kHz)",
                        ylabel="R1rho (s⁻¹)")
    hlines!(ax_fit_R1rho, [0]; linewidth=0)
    errorbars!(ax_fit_R1rho, state[:expfiterror]; color=c2)
    scatter!(ax_fit_R1rho, state[:expfitpoints]; label="Exponential fits", color=c2)
    lines!(ax_fit_R1rho, state[:fitR1rho]; label="Global fit")
    axislegend(ax_fit_R1rho; position=:rt)
    save(joinpath(outputdir, "dispersion.pdf"), fig; backend=CairoMakie)

    # intensities
    for i in 1:state[:nseries]
        title = "Peak integrals (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz)"
        filename = "intensities_$(round(0.001*νSL(state[:dataset])[i],digits=2))_kHz.pdf"
        fig = Figure()
        ax_fit = Axis(fig[1, 1];
                      xlabel="TSL (ms)",
                      ylabel="Peak integral",
                      title=title)
        errorbars!(ax_fit, state[:errorpoints][][i])
        scatter!(ax_fit, state[:scatterpoints][][i]; label="Observed")
        lines!(ax_fit, state[:fitseries][][i]; label="Global fit")
        axislegend(ax_fit; position=:rt)
        save(joinpath(outputdir, filename), fig; backend=CairoMakie)
    end

    # save fit results to CSV: peak/noise positions, and fit parameters. don't use additional libraries
    filename = joinpath(outputdir, "results.txt")
    open(filename, "w") do f
        println(f, "R1rho fitting results")
        println(f, "=====================")
        println(f, "")
        for filename in state[:filenames]
            println(f, "Input file: $filename")
        end
        println(f, "")
        println(f, "Peak position (ppm): $(state[:peakppm][])")
        println(f, "Noise position (ppm): $(state[:noiseppm][])")
        println(f, "Integration width (ppm): $(state[:dx][])")
        println(f, "")
        # report if exchange is significant
        f_stat, p_value = state[:ftest][]
        if p_value < 0.05
            println(f, "Exchange is significant (p < 0.05)")
        else
            println(f, "Exchange is not significant (p ≥ 0.05)")
        end
        println(f, "F-statistic: $(round(f_stat,digits=2))")
        println(f, "p-value: $(round(p_value,digits=4))")
        println(f, "")
        println(f, "Initial I0: $(state[:initialI0][])")
        println(f, "Initial R2,0 (s⁻¹): $(state[:initialR20][])")
        println(f, "Initial Rex (s⁻¹): $(state[:initialRex][])")
        println(f, "Initial kex (s⁻¹): $(exp(state[:initiallnk][]))")
        println(f, "Δδ stdev (ppm): $(state[:σΔδ][])")
        println(f, "")
        println(f, "Fitted I0: $(state[:fitI0][])")
        println(f, "Fitted R2,0 (s⁻¹): $(state[:fitR20][])")
        println(f, "Fitted Rex (s⁻¹): $(state[:fitRex][])")
        println(f, "Fitted K (s⁻¹): $(state[:fitK][])")
        println(f, "Fitted kex (s⁻¹): $(state[:fitkex][])")
        # include no-exchange fit result
        println(f, "")
        println(f, "Fitted R2,0 no-exchange (s⁻¹): $(state[:fitR20_null][])")
    end

    # write dispersion fit data to CSVs
    filename = joinpath(outputdir, "dispersion-fit.csv")
    open(filename, "w") do f
        println(f, "vSL (kHz),R1rho fit (s-1)")
        for i in state[:fitR1rho][]
            println(f, "$(i[1]),$(i[2])")
        end
    end
    filename = joinpath(outputdir, "dispersion-points.csv")
    open(filename, "w") do f
        println(f, "vSL (kHz),R1rho (s-1) [exponential fit],error (s-1)")
        for i in state[:expfiterror][]
            println(f, "$(i[1]),$(i[2]),$(i[3])")
        end
    end

    # write intensity data and fits to CSVs
    for i in 1:state[:nseries]
        filename = "intensities_$(round(0.001*νSL(state[:dataset])[i],digits=2))_kHz-fit.csv"
        filename = joinpath(outputdir, filename)
        open(filename, "w") do f
            println(f, "TSL (ms),fitted intensity")
            for p in state[:fitseries][][i]
                println(f, "$(p[1]),$(p[2])")
            end
        end
        filename = "intensities_$(round(0.001*νSL(state[:dataset])[i],digits=2))_kHz-points.csv"
        filename = joinpath(outputdir, filename)
        open(filename, "w") do f
            println(f, "TSL (ms),intensity,error")
            for p in state[:errorpoints][][i]
                println(f, "$(p[1]),$(p[2]),$(p[3])")
            end
        end
    end
end
