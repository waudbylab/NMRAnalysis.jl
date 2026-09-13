"Update experiment with non-overlapping clusters for fitting (e.g. [[1], [2, 3]])"
function cluster!(expt)
    adjacency = makeadjacency(expt.peaks[], expt)
    expt.clusters[] = connected_components(SimpleGraph(adjacency))
    @debug "$(length(expt.clusters[])) clusters of peaks found" expt.clusters[] maxlog = 10
end

"Generate adjacency matrix for peak overlap"
function makeadjacency(peaks, expt)
    npeaks = length(peaks)
    adjacency = zeros(npeaks, npeaks)

    for i in 1:npeaks
        for j in (i + 1):npeaks
            adjacency[i, j] = isoverlapping(peaks[i], peaks[j], expt)
            adjacency[j, i] = adjacency[i, j]
        end
    end

    return adjacency
end
