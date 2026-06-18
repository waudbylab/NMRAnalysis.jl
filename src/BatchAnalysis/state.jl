function initialisestate(spectra; referenceshift=0.0, standardx=0.0, standarddx=0.05,
                         unknownx=7.0, unknowndx=0.05)
    state = Dict{String,Any}()

    # spectra
    state["spectra"] = Observable(spectra)
    state["nspectra"] = lift(length, state["spectra"])
    # state["labels"] = lift(label, state["spectra"])
    state["selected"] = Observable(fill(true, state["nspectra"][]))
    state["nselected"] = lift(sum, state["selected"])
    state["selectedspectra"] = lift((spectra, selected) -> spectra[selected],
                                    state["spectra"], state["selected"])

    # chemical shift referencing
    state["referenced"] = Observable(false)
    state["referenceshift"] = Observable(referenceshift)

    # integration
    state["integrating"] = Observable(false)
    state["standardx"] = Observable(standardx)
    state["standarddx"] = Observable(standarddx)
    state["unknownx"] = Observable(unknownx)
    state["unknowndx"] = Observable(unknowndx)

    state["standardx1"] = lift((x, dx) -> x - dx / 2, state["standardx"],
                               state["standarddx"])
    state["standardx2"] = lift((x, dx) -> x + dx / 2, state["standardx"],
                               state["standarddx"])
    state["unknownx1"] = lift((x, dx) -> x - dx / 2, state["unknownx"], state["unknowndx"])
    state["unknownx2"] = lift((x, dx) -> x + dx / 2, state["unknownx"], state["unknowndx"])

    state["standardintegrals"] = lift(integrate, state["spectra"], state["standardx1"],
                                      state["standardx2"])
    state["unknownintegrals"] = lift(integrate, state["spectra"], state["unknownx1"],
                                     state["unknownx2"])
    state["integralratios"] = lift((x, ref) -> x ./ ref, state["unknownintegrals"],
                                   state["standardintegrals"])

    # filenames
    state["figurefilename"] = Observable("spectra.pdf")
    state["resultsfilename"] = Observable("results.csv")

    return state
end

function preparegui!(state)
    GLMakie.activate!()

    state["gui"] = Dict{String,Any}()
    gui = state["gui"]

    # spectrum plotting
    gui["plotscale"] = Observable(1.0)
    gui["linedata"] = lift(linedata, state["selectedspectra"], gui["plotscale"])
    gui["linecolors"] = lift(linecolors, state["nselected"])

    # integration
    gui["integration_label_1"] = lift((x, dx) -> "Standard: $(round(x,digits=3)) ± $(round(dx,digits=3)) ppm",
                                      state["standardx"], state["standarddx"])
    gui["integration_label_2"] = lift((x, dx) -> "Unknown: $(round(x,digits=3)) ± $(round(dx,digits=3)) ppm",
                                      state["unknownx"], state["unknowndx"])
    gui["integralplate"] = lift(ratios -> reshape(ratios, 12, 8), state["integralratios"])

    # filenames
    gui["spectra_filename"] = Observable("spectra-plot.pdf")
    gui["results_filename"] = Observable("results.csv")

    # selected points for heatmap
    return gui["selected_points"] = lift(selectedpoints, state["selected"])
end

function selectedpoints(selected)
    p = [Point2f(i, j) for i in 1:12, j in 1:8]
    return p[selected]
end