"""
    STDExperiment(dataset; regions, reference=:reference, excess=1.0)

Saturation-transfer difference analysis. The planes carry two arrayed variables:

- `sat`  : the saturation condition — one designated `reference` (off-resonance) value
           plus one or more on-resonance saturation frequencies
           (e.g. `:reference, :methyl, :aromatic`).
- `tsat` : the saturation time (s).

For each ligand `region`, each non-reference saturation `sat`, and each `tsat`, the STD
fraction is contrasted against the reference at matching `tsat`:

    STD = (I_reference − I_sat) / I_reference

and the STD amplification factor `STD-AF = STD · excess`, where `excess = [L]/[P]`
(default 1, i.e. report STD%). When several `tsat` values are present for a
(region, sat) pair, a buildup curve `STD-AF(t) = STD-AF_max·(1 − exp(−k·t))` is fitted
and the initial slope `STD-AF₀ = STD-AF_max·k` reported (this removes the T1 bias that
makes raw STD% an unreliable epitope ranking). Finally an epitope map normalises STD-AF
across regions to the strongest signal.

This is a *contrast* experiment (the 1D analogue of GUI2D's `hetnoe2d`); it overrides
`analyse` rather than using the curve-fit pipeline.
"""
struct STDExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    reference::Any
    excess::Float64
end

function STDExperiment(dataset::Dataset1D; regions, reference=:reference, excess::Real=1.0)
    haskey(first(dataset.planes.vars), :sat) ||
        throw(ArgumentError("STD planes must carry a :sat variable"))
    haskey(first(dataset.planes.vars), :tsat) ||
        throw(ArgumentError("STD planes must carry a :tsat variable"))
    reference in column(dataset.planes, :sat) ||
        throw(ArgumentError("reference saturation $(repr(reference)) not present in :sat"))
    return STDExperiment(dataset, collect(Region, regions), reference, Float64(excess))
end

dataset(e::STDExperiment) = e.dataset
regions(e::STDExperiment) = e.regions

"""STD fraction (× excess) for one region, saturation and saturation time."""
struct STDPoint
    region::String
    sat::Any
    tsat::Float64
    std::Measurement{Float64}
end

"""Buildup fit for one region/saturation: initial slope, plateau and rate."""
struct STDBuildup
    region::String
    sat::Any
    std_af0::Measurement{Float64}
    std_af_max::Measurement{Float64}
    k::Measurement{Float64}
    converged::Bool
end

"""Epitope entry: STD-AF and its value relative to the strongest region."""
struct EpitopePoint
    region::String
    sat::Any
    tsat::Float64
    std_af::Measurement{Float64}
    relative::Float64
end

function BuildupModel()
    return CurveFitModel((x, p) -> @.(p[1] * (1 - exp(-p[2] * x))),
                         ["STD_AF_max", "k"],
                         (x, y) -> [maximum(y), 3.0 / maximum(x)];
                         xlabel="Saturation time / s", ylabel="STD-AF")
end

"""
    analyse(e::STDExperiment, [dataset, regions]) -> NamedTuple

Returns `(; points, buildups, epitope)`. As for the curve-fit experiments, the dataset
and regions may be supplied explicitly so the GUI can recompute live.
"""
analyse(e::STDExperiment) = analyse(e, e.dataset, e.regions)

function analyse(e::STDExperiment, ds::Dataset1D, regs)
    sats = unique(column(ds.planes, :sat))
    onres = filter(!=(e.reference), sats)

    # Per-region integrals over every plane.
    points = STDPoint[]
    buildups = STDBuildup[]
    for region in regs
        I = integrals(region, ds)
        sat_of = column(ds.planes, :sat)
        tsat_of = column(ds.planes, :tsat)

        # Reference integral at each saturation time (average duplicates).
        reftimes = unique(tsat_of[sat_of .== e.reference])
        refI = Dict(τ => mean(I[(sat_of .== e.reference) .& (tsat_of .== τ)])
                    for τ in reftimes)

        for s in onres
            stds = STDPoint[]
            for τ in sort(unique(tsat_of[sat_of .== s]))
                haskey(refI, τ) || continue
                Is = mean(I[(sat_of .== s) .& (tsat_of .== τ)])
                std = e.excess * (refI[τ] - Is) / refI[τ]
                pt = STDPoint(region.label, s, τ, std)
                push!(points, pt)
                push!(stds, pt)
            end
            # Buildup fit when enough saturation times are available.
            if length(stds) ≥ 3
                x = [p.tsat for p in stds]
                y = [p.std for p in stds]
                fit = fit_series(BuildupModel(), x, y)
                stdmax = fit.params[1]
                k = fit.params[2]
                push!(buildups,
                      STDBuildup(region.label, s, stdmax * k, stdmax, k, fit.converged))
            end
        end
    end

    epitope = epitope_map(points)
    return (; points, buildups, epitope)
end

"""
    epitope_map(points) -> Vector{EpitopePoint}

Normalise STD-AF across regions, separately for each (sat, tsat), to the strongest
region.
"""
function epitope_map(points::Vector{STDPoint})
    epitope = EpitopePoint[]
    for s in unique(p.sat for p in points)
        for τ in unique(p.tsat for p in points if p.sat == s)
            group = filter(p -> p.sat == s && p.tsat == τ, points)
            maxstd = maximum(Measurements.value(p.std) for p in group)
            maxstd == 0 && continue
            for p in group
                rel = Measurements.value(p.std) / maxstd
                push!(epitope, EpitopePoint(p.region, s, τ, p.std, rel))
            end
        end
    end
    return epitope
end
