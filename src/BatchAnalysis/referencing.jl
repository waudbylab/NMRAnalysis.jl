function referencespectra!(state, xl1, xl2)
    refshift = state["referenceshift"][]

    state["spectra"][] = map(state["spectra"][]) do spec
        xs = data(spec[xl1 .. xl2], F1Dim)
        ys = data(spec[xl1 .. xl2])
        # find position of maximum intensity
        _, i = findmax(ys)
        peakposition = xs[i]
        # shift spectrum to move peak to reference shift
        # dx = refshift - peakposition
        # add_offset(spec, F1Dim, dx)
        return spec = reference(spec, F1Dim, peakposition => refshift)
    end
end