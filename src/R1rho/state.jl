function initialisestate(dataset)
    state = Dict{Symbol,Any}()
    state[:dataset] = dataset
    state[:onres] = isonres(dataset)

    # Find the maximum of dataset.spectra[1]
    max_idx = argmax(dataset.spectra[1])
    max_chemical_shift = data(dataset.spectra[1], F1Dim)[max_idx]

    # Estimate a noise position from the point 90% of the way across the spectrum
    noise_idx = round(Int, 0.9 * length(dataset.spectra[1]))
    noise_chemical_shift = data(dataset.spectra[1], F1Dim)[noise_idx]

    state[:peakppm] = Observable(max_chemical_shift)
    state[:noiseppm] = Observable(noise_chemical_shift)
    state[:dx] = Observable(0.5) # ppm
    state[:σΔδ] = Observable(2.0) # ppm

    state[:intensities] = lift((x, dx) -> integrate(dataset, x, dx), state[:peakppm],
                               state[:dx])
    state[:noise] = lift((x, dx) -> noise(dataset, x, dx), state[:noiseppm], state[:dx])

    state[:isfitting] = Observable(true)
    state[:outputdir] = Observable("out")

    if state[:onres]
        # data points
        state[:series] = onresseries(dataset)
        state[:nseries] = length(state[:series])
        state[:scatterpoints] = lift(state[:intensities]) do I
            return [Point2f.(1000 * dataset.TSLs[state[:series][i]], I[state[:series][i]])
                    for i in 1:state[:nseries]]
        end
        state[:errorpoints] = lift(state[:intensities], state[:noise]) do I, σ
            return [tuple.(1000 * dataset.TSLs[state[:series][i]], I[state[:series][i]], σ)
                    for i in 1:state[:nseries]]
        end

        # current series
        state[:currentseries] = Observable(1)
        state[:currentscatter] = lift(state[:currentseries],
                                      state[:scatterpoints]) do i, points
            return points[i]
        end
        state[:currenterror] = lift(state[:currentseries], state[:errorpoints]) do i, points
            return points[i]
        end
        state[:currentspectrum] = lift(state[:currentseries]) do i
            j = state[:series][i][1]
            return Point2f.(data(dataset.spectra[j], F1Dim), data(dataset.spectra[j]))
        end

        # parameters
        state[:initialI0] = Observable(state[:intensities][][1])
        # fit exponentials to all series to get range of initial R2 values - use min and max to set initial R20 and Rex
        quick_R2_vals = pmean.([fitexp(dataset.TSLs[state[:series][i]],
                                       state[:intensities][][state[:series][i]],
                                       [state[:initialI0][]])
                                for i in 1:state[:nseries]])
        quick_initial_R20 = minimum(quick_R2_vals)
        if quick_initial_R20 < 1.0
            quick_initial_R20 = 1.0
        end
        if quick_initial_R20 > 100.0
            quick_initial_R20 = 100.0
        end
        quick_initial_Rex = maximum(quick_R2_vals) - quick_initial_R20
        if quick_initial_Rex > 100.0
            quick_initial_Rex = 100.0
        end
        if quick_initial_Rex < 1.1
            quick_initial_Rex = 0.1
        end
        state[:initialR20] = Observable(quick_initial_R20)
        state[:initialRex] = Observable(quick_initial_Rex)
        # set initial k to π * maximum vSL
        max_vSL = maximum(νSL(dataset))
        state[:initiallnk] = Observable(log(π * max_vSL))
        state[:initialpars] = lift(state[:initialI0], state[:initialR20],
                                   state[:initialRex],
                                   state[:initiallnk]) do I0, R20, Rex, lnk
            return [I0, R20, Rex, lnk]
        end

        # fitting
        state[:fit] = lift(state[:isfitting], state[:initialpars],
                           state[:intensities]) do isfitting, p0, _
            # if isfitting
            return fit_onres(state, p0)
            # else
            #     LsqFit.LsqFitResult([1.0, 2.0, 10.0, 8.5],Float64[],zeros(4,4),false,Vector{LsqFit.LMState{LsqFit.LevenbergMarquardt}}(),Float64[])
            # end
        end
        # fit results
        state[:fitpars] = lift(coef, state[:fit])
        state[:fiterrs] = lift(state[:fit]) do fit
            try
                stderrors(fit)
            catch
                [0.1, 0.1, 0.1, 0.1]
            end
        end
        state[:fitI0] = lift((p, e) -> p[1] ± e[1], state[:fitpars], state[:fiterrs])
        state[:fitR20] = lift((p, e) -> p[2] ± e[2], state[:fitpars], state[:fiterrs])
        state[:fitRex] = lift((p, e) -> p[3] ± e[3], state[:fitpars], state[:fiterrs])
        state[:fitlnk] = lift((p, e) -> p[4] ± e[4], state[:fitpars], state[:fiterrs])
        state[:fitR2s] = lift(state[:fitpars]) do p
            return [fitexp(dataset.TSLs[state[:series][i]],
                           state[:intensities][][state[:series][i]], p)
                    for i in 1:state[:nseries]]
        end

        # extract K and kex
        bf = dataset.spectra[1][1, :bf] # larmor frequency
        state[:fitK] = lift(lnK -> exp(lnK), state[:fitlnk])
        state[:fitkex] = lift((K, σΔδ) -> estimatekexfromK(K, σΔδ, bf),
                              state[:fitK], state[:σΔδ])

        # null fit
        state[:fit_null] = lift(state[:isfitting], state[:initialpars],
                                state[:intensities]) do isfitting, p0, _
            return fit_onres_null(state, p0[1:2])
        end
        state[:fitpars_null] = lift(coef, state[:fit_null])
        state[:fiterrs_null] = lift(state[:fit_null]) do fit
            try
                stderrors(fit)
            catch
                [0.1, 0.1]
            end
        end
        state[:fitI0_null] = lift((p, e) -> p[1] ± e[1], state[:fitpars_null],
                                  state[:fiterrs_null])
        state[:fitR20_null] = lift((p, e) -> p[2] ± e[2], state[:fitpars_null],
                                   state[:fiterrs_null])

        # Model comparison statistics
        state[:ftest] = lift(state[:fit], state[:fit_null],
                             state[:intensities]) do full_fit, null_fit, intensities
            n = length(intensities)
            p_full = 4  # Parameters in full model
            p_null = 2  # Parameters in null model

            rss_full = sum(abs2, full_fit.resid)
            rss_null = sum(abs2, null_fit.resid)

            f_stat = ((rss_null - rss_full) / (p_full - p_null)) / (rss_full / (n - p_full))
            p_value = 1 - cdf(FDist(p_full - p_null, n - p_full), f_stat)

            return (f_stat, p_value)
        end

        # fit plotting
        state[:fitseries] = lift(state[:fitpars]) do p
            ts = range(0, 1.1 * maximum(TSL(dataset)), 100)
            return [Point2f.(1000 * ts, model_I_onres(ts, νSL(dataset)[i], p))
                    for i in 1:state[:nseries]]
        end
        state[:currentfit] = lift(state[:currentseries], state[:fitseries]) do i, points
            return points[i]
        end
        state[:fitR1rho] = lift(state[:fitpars]) do p
            xs = range(0, 1.1 * maximum(νSL(dataset)), 100)
            return Point2f.(0.001 * xs, model_R1rho_onres(xs, p))
        end
        state[:expfitpoints] = lift(state[:fitR2s]) do R2s
            return Point2f.(0.001 * νSL(dataset), pmean.(R2s))
        end
        state[:expfiterror] = lift(state[:fitR2s]) do R2s
            return tuple.(0.001 * νSL(dataset), pmean.(R2s), pstd.(R2s))
        end
        state[:currentexpfit] = lift(state[:currentseries], state[:fitR2s]) do i, R2s
            ts = range(0, 1.1 * maximum(TSL(dataset)), 100)
            R2 = pmean.(R2s)[i]
            I0 = state[:fitpars][][1]
            ys = @. I0 * exp(-ts * R2)
            return Point2f.(1000 * ts, ys)
        end
        state[:residualpoints] = lift(state[:currentseries], state[:intensities],
                                      state[:fitpars]) do i, I, p
            t = dataset.TSLs[state[:series][i]]
            y = I[state[:series][i]]
            yfit = model_I_onres(t, νSL(dataset)[i], p)
            return Point2f.(1000 * t, y - yfit)
        end
        state[:residualerror] = lift(state[:currentseries], state[:intensities],
                                     state[:noise], state[:fitpars]) do i, I, σ, p
            t = dataset.TSLs[state[:series][i]]
            y = I[state[:series][i]]
            ye = σ
            yfit = model_I_onres(t, νSL(dataset)[i], p)
            return tuple.(1000 * t, y - yfit, ye)
        end

        # null fit plotting
        state[:fitseries_null] = lift(state[:fitpars_null]) do p
            ts = range(0, 1.1 * maximum(TSL(dataset)), 100)
            return [Point2f.(1000 * ts, model_I_onres_null(ts, νSL(dataset)[i], p))
                    for i in 1:state[:nseries]]
        end
        state[:currentfit_null] = lift(state[:currentseries],
                                       state[:fitseries_null]) do i, points
            return points[i]
        end
        state[:fitR1rho_null] = lift(state[:fitpars_null]) do p
            xs = range(0, 1.1 * maximum(νSL(dataset)), 100)
            return Point2f.(0.001 * xs, model_R1rho_onres_null(xs, p))
        end
        state[:residualpoints_null] = lift(state[:currentseries], state[:intensities],
                                           state[:fitpars_null]) do i, I, p
            t = dataset.TSLs[state[:series][i]]
            y = I[state[:series][i]]
            yfit = model_I_onres_null(t, νSL(dataset)[i], p)
            return Point2f.(1000 * t, y - yfit)
        end
        state[:residualerror_null] = lift(state[:currentseries], state[:intensities],
                                          state[:noise], state[:fitpars_null]) do i, I, σ, p
            t = dataset.TSLs[state[:series][i]]
            y = I[state[:series][i]]
            ye = σ
            yfit = model_I_onres_null(t, νSL(dataset)[i], p)
            return tuple.(1000 * t, y - yfit, ye)
        end
    end

    return state
end
