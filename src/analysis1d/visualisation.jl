# Result-panel visualisation. Pure data builders (result -> plot primitives) dispatched on
# the experiment type, plus axis labels and a summary string. The GUI lifts Observables off
# these; saving reuses the same builders, so live and exported plots share one code path.

"""
    result_plotdata(expt, result, activelabel) -> (points, errors, fitlines)

Plot primitives for the active region: observed `points` (`Vector{Point2f}`), `errors`
(`Vector{NTuple{3,Float64}}` of `(x, y, σ)`), and `fitlines` (`Vector{Point2f}`, series
separated by `NaN`). The generic method covers the curve-fit / NoFitting experiments;
STD overrides it.
"""
function result_plotdata(::Experiment1D, result, activelabel::AbstractString)
    series = filter(s -> s.region == activelabel, result.series)
    points = Point2f[]
    errors = Tuple{Float64,Float64,Float64}[]
    fitlines = Point2f[]
    for s in series
        for k in eachindex(s.x)
            v = Measurements.value(s.y[k])
            e = Measurements.uncertainty(s.y[k])
            push!(points, Point2f(s.x[k], v))
            push!(errors, (s.x[k], v, e))
        end
        if !(s.model isa NoFitting) && !isempty(s.params)
            xs = collect(range(min(0.0, minimum(s.x)), 1.05 * maximum(s.x), 100))
            ys = s.model.func(xs, Measurements.value.(s.params))
            append!(fitlines, Point2f.(xs, ys))
            push!(fitlines, Point2f(NaN32, NaN32))
        end
    end
    return (points, errors, fitlines)
end

function result_plotdata(::STDExperiment, result, activelabel::AbstractString)
    pts = filter(p -> p.region == activelabel, result.points)
    points = [Point2f(p.tsat, Measurements.value(p.std)) for p in pts]
    errors = [(p.tsat, Measurements.value(p.std), Measurements.uncertainty(p.std))
              for p in pts]
    fitlines = Point2f[]
    tmax = isempty(pts) ? 1.0 : maximum(p.tsat for p in pts)
    for b in filter(b -> b.region == activelabel, result.buildups)
        xs = collect(range(0.0, 1.05 * tmax, 100))
        smax = Measurements.value(b.std_af_max)
        k = Measurements.value(b.k)
        ys = @. smax * (1 - exp(-k * xs))
        append!(fitlines, Point2f.(xs, ys))
        push!(fitlines, Point2f(NaN32, NaN32))
    end
    return (points, errors, fitlines)
end

# axis labels for the result panel
result_labels(::Experiment1D) = ("x", "Integral")
result_labels(::RelaxationExperiment) = ("Relaxation delay / s", "Integral")
result_labels(::TractExperiment) = ("Relaxation delay / s", "Integral")
result_labels(::NutationExperiment) = ("Pulse duration / s", "Integral")
result_labels(::KineticsExperiment) = ("Time", "Integral")
result_labels(::STDExperiment) = ("Saturation time / s", "STD fraction")

# ---- summary text -------------------------------------------------------------

"""
    summary_text(expt, result) -> String

Human-readable summary of the fit, shown in the GUI and written to `summary.txt`.
"""
function summary_text(e::Experiment1D, result)
    io = IOBuffer()
    for s in result.series
        grp = isempty(s.group) ? "" : " " * string(s.group)
        println(io, "$(s.region)$grp:")
        for (name, p) in zip(s.names, s.params)
            println(io, "  $name = $p")
        end
    end
    _summary_extra(io, e, result.summary)
    return String(take!(io))
end

_summary_extra(::IOBuffer, ::Experiment1D, ::Nothing) = nothing
function _summary_extra(io::IOBuffer, ::TractExperiment, summary)
    println(io, "")
    for s in summary
        println(io, "$(s.region): τc = $(s.τc) ns   (ΔR = $(s.Ranti - s.Rtrosy) s⁻¹)")
    end
end
function _summary_extra(io::IOBuffer, ::NutationExperiment, summary)
    println(io, "")
    for s in summary
        println(io, "$(s.region): ν₁ = $(s.ν) Hz, 90° = $(1e6 * s.pulse90) µs")
    end
end

function summary_text(::STDExperiment, result)
    io = IOBuffer()
    println(io, "STD results:")
    for p in result.points
        println(io, "  $(p.region) / $(p.sat) @ $(p.tsat) s: STD = $(p.std)")
    end
    if !isempty(result.buildups)
        println(io, "\nBuildup (initial slope):")
        for b in result.buildups
            println(io, "  $(b.region) / $(b.sat): STD-AF₀ = $(b.std_af0), k = $(b.k) s⁻¹")
        end
    end
    if !isempty(result.epitope)
        println(io, "\nEpitope (relative):")
        for ep in result.epitope
            println(io, "  $(ep.region) / $(ep.sat): $(round(100 * ep.relative; digits=0)) %")
        end
    end
    return String(take!(io))
end
