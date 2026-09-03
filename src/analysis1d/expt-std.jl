# STD: saturation-transfer difference, with buildup curves and an epitope map.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    std1d(spec, sat, tsat; regions=nothing, reference=:reference, excess=1.0, integration=nothing)

Analyse an STD experiment. `sat` and `tsat` are per-plane vectors giving the saturation
condition (one of which is the `reference`) and the saturation time. Reports STD fractions
per region, buildup curves where several saturation times are present, and an epitope map.

`regions` defaults to a single peak-detected region; add further named ligand regions
interactively in the GUI (press `A`) rather than specifying them up front.
"""
function std1d(spec, sat::AbstractVector, tsat::AbstractVector; regions=nothing,
               reference=:reference, excess::Real=1.0, integration=nothing)
    spec = _spec(spec)
    vars = [(; sat=sat[i], tsat=Float64(tsat[i])) for i in eachindex(sat)]
    ds = dataset_from_spec(spec, vars)
    expt = isnothing(regions) ? STDExperiment(ds; reference, excess) :
           STDExperiment(ds; regions, reference, excess)
    return run1d(expt; integration)
end

# ---- 2. type ------------------------------------------------------------------

"""
    STDExperiment(dataset; regions=[defaultregion(dataset)], reference=:reference, excess=1.0)

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

function STDExperiment(dataset::Dataset1D; regions=[defaultregion(dataset)],
                       reference=:reference, excess::Real=1.0)
    haskey(first(dataset.planes.vars), :sat) ||
        throw(ArgumentError("STD planes must carry a :sat variable"))
    haskey(first(dataset.planes.vars), :tsat) ||
        throw(ArgumentError("STD planes must carry a :tsat variable"))
    reference in column(dataset.planes, :sat) ||
        throw(ArgumentError("reference saturation $(repr(reference)) not present in :sat"))
    return STDExperiment(dataset, collect(Region, regions), reference, Float64(excess))
end

# ---- 3. interface -------------------------------------------------------------
# STD overrides `analyse` wholesale, so it declares no fitaxis/groupcols.

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

# ---- 4. science ---------------------------------------------------------------

"""
    ContrastModel(; reference)

A categorical series model: within a group of slices it contrasts each non-reference
slice against the `reference` slice. The contrast value is `(I_ref − I_slice)/I_ref`
(the STD fraction). Used by the STD experiment; the 1D analogue of GUI2D's `hetnoe2d`.
"""
struct ContrastModel <: SeriesModel
    reference::Any
end

ContrastModel(; reference) = ContrastModel(reference)

function BuildupModel()
    return CurveFitModel((x, p) -> @.(p[1] * (1 - exp(-p[2] * x))),
                         ["STD_AF_max", "k"],
                         (x, y) -> [maximum(y), 3.0 / maximum(x)];
                         xlabel="Saturation time / s", ylabel="STD-AF")
end

"""
    analyse(e::STDExperiment, [dataset, regions]; isfitting=true) -> NamedTuple

Returns `(; points, buildups, epitope)`. As for the curve-fit experiments, the dataset
and regions may be supplied explicitly so the GUI can recompute live. `isfitting=false`
(the GUI's Fitting toggle switched off) skips the buildup curve fit - the raw STD
fractions (`points`) still come through, but `buildups` (and, following from it,
`epitope`) come back empty.
"""
analyse(e::STDExperiment) = analyse(e, e.dataset, e.regions)

function analyse(e::STDExperiment, ds::Dataset1D, regs; isfitting::Bool=true)
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
            if isfitting && length(stds) ≥ 3
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

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::STDExperiment) = "STD"

resultlabels(::STDExperiment) = ("Saturation time / s", "STD fraction")

spectruminfo(::STDExperiment, vars::NamedTuple) = "$(vars.sat), $(round(vars.tsat; digits=3)) s sat"

function resultplotdata(::STDExperiment, result, activelabel::AbstractString)
    pts = filter(p -> p.region == activelabel, result.points)
    sats = unique(p.sat for p in pts)
    return map(sats) do sat
        ps = filter(p -> p.sat == sat, pts)
        points = [Point2f(p.tsat, Measurements.value(p.std)) for p in ps]
        errors = [(p.tsat, Measurements.value(p.std), Measurements.uncertainty(p.std))
                  for p in ps]
        bi = findfirst(b -> b.region == activelabel && b.sat == sat, result.buildups)
        fitline = if !isnothing(bi)
            b = result.buildups[bi]
            tmax = maximum(p.tsat for p in ps)
            xs = collect(range(0.0, 1.05 * tmax, 100))
            smax, k = Measurements.value(b.std_af_max), Measurements.value(b.k)
            Point2f.(xs, @.(smax * (1 - exp(-k * xs))))
        else
            Point2f[]
        end
        ResultSeries(points, errors, fitline, string(sat))
    end
end

# STD's result has no `.series`/`.summary` fields to split into a resultsheader/
# secondarytext pair (it's a contrast rather than a curve-fit experiment) - show its
# whole summary in the secondary block instead.
resultsheader(::STDExperiment, result, activelabel::AbstractString) = ""
function secondarytext(e::STDExperiment, result, activelabel::AbstractString)
    return summary_text(e, result, activelabel)
end

function summary_text(::STDExperiment, result, activelabel::AbstractString)
    io = IOBuffer()
    println(io, "STD fractions")
    println(io, "-------------")
    for p in result.points
        p.region == activelabel || continue
        println(io, "  $(p.region) / $(p.sat) @ $(p.tsat) s : $(fmt(p.std))")
    end
    buildups = filter(b -> b.region == activelabel, result.buildups)
    if !isempty(buildups)
        println(io, "\nBuildup (initial slope, T1-corrected)")
        println(io, "--------------------------------------")
        for b in buildups
            println(io,
                    "  $(b.region) / $(b.sat) : STD-AF₀ = $(fmt(b.std_af0)), k = $(fmt(b.k)) s⁻¹")
        end
    end
    epitope = filter(ep -> ep.region == activelabel, result.epitope)
    if !isempty(epitope)
        println(io, "\nEpitope map (relative to strongest signal)")
        println(io, "-------------------------------------------")
        for ep in epitope
            println(io, "  $(ep.region) / $(ep.sat) : $(round(100 * ep.relative; digits=0)) %")
        end
    end
    return String(take!(io))
end
