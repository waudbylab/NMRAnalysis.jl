# Result-panel visualisation. Pure data builders (result -> plot primitives) dispatched on
# the experiment type, plus axis labels and a summary string. The GUI lifts Observables off
# these; saving reuses the same builders, so live and exported plots share one code path.

"""
    ResultSeries

One coloured curve in the result panel: observed `points` (`Vector{Point2f}`), `errors`
(`Vector{NTuple{3,Float64}}` of `(x, y, σ)`), a `fitline` (`Vector{Point2f}`, empty if
unfitted), and a `label` (e.g. the TRACT `which`, or the STD saturation frequency; empty
for an ungrouped single series).
"""
struct ResultSeries
    points::Vector{Point2f}
    errors::Vector{Tuple{Float64,Float64,Float64}}
    fitline::Vector{Point2f}
    label::String
end

"""
    result_plotdata(expt, result, activelabel) -> Vector{ResultSeries}

Plot primitives for the active region, one `ResultSeries` per group (e.g. TRACT's
TROSY/anti-TROSY pair, or STD's saturation frequencies) so each can be drawn in a
distinct, matching colour. The generic method covers the curve-fit / NoFitting
experiments; STD overrides it.
"""
function result_plotdata(::Experiment1D, result, activelabel::AbstractString)
    series = filter(s -> s.region == activelabel, result.series)
    return map(series) do s
        points = Point2f.(s.x, Measurements.value.(s.y))
        errors = [(s.x[k], Measurements.value(s.y[k]), Measurements.uncertainty(s.y[k]))
                  for k in eachindex(s.x)]
        fitline = if !(s.model isa NoFitting) && !isempty(s.params)
            xs = collect(range(min(0.0, minimum(s.x)), 1.05 * maximum(s.x), 100))
            Point2f.(xs, s.model.func(xs, Measurements.value.(s.params)))
        else
            Point2f[]
        end
        ResultSeries(points, errors, fitline, groupname(s.group))
    end
end

function result_plotdata(::STDExperiment, result, activelabel::AbstractString)
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

"""Human-readable label for a grouping key, e.g. `(which = :trosy,)` → `"trosy"`."""
function groupname(group::NamedTuple)
    isempty(group) && return ""
    return join((string(v) for v in values(group)), ", ")
end

# axis labels for the result panel
result_labels(::Experiment1D) = ("x", "Intensity / noise")
result_labels(::RelaxationExperiment) = ("Relaxation delay / s", "Intensity / noise")
result_labels(::TractExperiment) = ("Relaxation delay / s", "Intensity / noise")
result_labels(::NutationExperiment) = ("Pulse duration / s", "Intensity / noise")
result_labels(::KineticsExperiment) = ("Time", "Intensity / noise")
result_labels(::DiffusionExperiment) = ("Relative gradient strength", "Intensity / noise")
result_labels(::STDExperiment) = ("Saturation time / s", "STD fraction")

# ---- summary text -------------------------------------------------------------

"""Compact, consistently-rounded rendering of a `Measurement` (or plain number)."""
fmt(x::Measurement, digits=4) = string(round(Measurements.value(x); sigdigits=digits),
                                       " ± ",
                                       round(Measurements.uncertainty(x); sigdigits=2))
fmt(x::Real, digits=4) = string(round(x; sigdigits=digits))
fmt(::Nothing, digits=4) = "n/a"

const PARAM_UNITS = Dict("R" => " s⁻¹", "ν" => " Hz")

"""
    summary_text(expt, result) -> String

Human-readable summary of the fit, shown in the GUI and written to `summary.txt`.

# Note
Output-format consistency across the different 1D/2D analyses (units, significant
figures, CSV vs text) is an open question — see `PLAN.md` — this only tidies the
single-experiment summary shown here.
"""
function summary_text(e::Experiment1D, result)
    io = IOBuffer()
    for s in result.series
        header = isempty(s.group) ? s.region : "$(s.region) ($(groupname(s.group)))"
        println(io, header)
        println(io, "-"^length(header))
        for (name, p) in zip(s.names, s.params)
            println(io, "  $(rpad(name, 6)) = $(fmt(p))$(get(PARAM_UNITS, name, ""))")
        end
        println(io)
    end
    _summary_extra(io, e, result.summary)
    return String(take!(io))
end

_summary_extra(::IOBuffer, ::Experiment1D, ::Nothing) = nothing
function _summary_extra(io::IOBuffer, ::TractExperiment, summary)
    for s in summary
        println(io, "$(s.region): τc = $(fmt(s.τc)) ns   (ΔR = $(fmt(s.Ranti - s.Rtrosy)) s⁻¹)")
    end
end
function _summary_extra(io::IOBuffer, ::NutationExperiment, summary)
    for s in summary
        println(io,
                "$(s.region): ν₁ = $(fmt(s.ν)) Hz, 90° pulse = $(fmt(1e6 * s.pulse90)) µs")
    end
end
function _summary_extra(io::IOBuffer, ::DiffusionExperiment, summary)
    for s in summary
        line = "$(s.region): D = $(fmt(s.D)) m² s⁻¹"
        isnothing(s.rH) || (line *= ", rH = $(fmt(s.rH)) Å (η = $(fmt(s.η)) mPa s)")
        println(io, line)
    end
end

function summary_text(::STDExperiment, result)
    io = IOBuffer()
    println(io, "STD fractions")
    println(io, "-------------")
    for p in result.points
        println(io, "  $(p.region) / $(p.sat) @ $(p.tsat) s : $(fmt(p.std))")
    end
    if !isempty(result.buildups)
        println(io, "\nBuildup (initial slope, T1-corrected)")
        println(io, "--------------------------------------")
        for b in result.buildups
            println(io,
                    "  $(b.region) / $(b.sat) : STD-AF₀ = $(fmt(b.std_af0)), k = $(fmt(b.k)) s⁻¹")
        end
    end
    if !isempty(result.epitope)
        println(io, "\nEpitope map (relative to strongest signal)")
        println(io, "-------------------------------------------")
        for ep in result.epitope
            println(io, "  $(ep.region) / $(ep.sat) : $(round(100 * ep.relative; digits=0)) %")
        end
    end
    return String(take!(io))
end
