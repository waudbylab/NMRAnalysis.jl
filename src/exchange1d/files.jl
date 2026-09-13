# Output files, written to the layout and column rules in
# `docs/src/advanced/conventions.md`:
#
#   out/
#     summary.txt        the record to read, including how to repeat the analysis
#     results.csv        one row per experiment: what was fitted, and its settings
#     series.csv         the measurements, one row per experiment per data point
#     global.csv         every fitted parameter, with its initial value and uncertainty
#     fit.pdf            all experiments on one grid
#     overlay_*.pdf      overlays of comparable experiments
#     experiments/<label>.pdf, experiments/<label>.csv
#
# The **entity** here is the experiment: a joint fit has no per-peak or per-region results,
# and every fitted parameter is shared between experiments (the model parameters globally,
# the spin parameters per field, the nuisance parameters per experiment type and field). So
# `results.csv` is an index of what went into the fit, and every number that came out of it
# is in `global.csv`. This is the same shape as a kinetics run, whose `results.csv` likewise
# carries keys and no fitted parameters.

# ---- per-experiment hooks -----------------------------------------------------

"""
    seriescoordinates(expt) -> Vector{Pair{Symbol,Any}}

The coordinates of each of `expt`'s data points, as `name => values`. A vector value gives
one entry per data point; a scalar is a setting constant across the experiment but varying
between experiments in the same problem (a saturation field, a spin-lock power), and is
repeated down the rows.

Coordinates are named per quantity rather than reduced to one `x` column, because a joint
fit mixes experiment types: a problem containing CEST and R1 experiments has both an
`offset (ppm)` and a `delay (s)` column, and each row fills in the ones that apply to it.
"""
seriescoordinates(e::R1Experiment) = [:delay => e.delays]
function seriescoordinates(e::CESTExperiment)
    return [:offset => e.δsat, :nu1 => e.ν1, :tsat => e.saturation_time]
end
seriescoordinates(e::R1rhoOnResExperiment) = [:nu_SL => e.νSL]
seriescoordinates(e::R1rhoOffResExperiment) = [:offset => e.offsets_ppm, :nu_SL => e.νSL]

"""
    observable(expt) -> (name, unit)

What `expt.observed_intensities` actually holds. CEST and R1 experiments observe an
integrated intensity, normalised to its own maximum and so dimensionless. Both R1ρ
experiments observe a *rate*, already fitted from the decay over the spin-lock delays, so
calling its column `I` would mislabel it.
"""
observable(::AbstractExperiment) = (:I, "")
observable(::R1rhoOnResExperiment) = (:R1rho, "s-1")
observable(::R1rhoOffResExperiment) = (:R1rho, "s-1")

const COORDINATE_UNITS = Dict(:delay => "s", :offset => "ppm", :nu1 => "Hz",
                              :tsat => "s", :nu_SL => "Hz")

coordinateunit(name::Symbol) = get(COORDINATE_UNITS, name, "")

"""
    parameterunit(item) -> String

ASCII unit for a fitted parameter, matched on its `ComponentArray` key with any `log`
prefix, section and state index already stripped by `_paramkey`.

Deliberately conservative: a key this does not recognise gets no unit rather than a guessed
one. `Kd` is the notable blank - it comes out in whatever units the sample metadata gives
concentrations in, which this module never learns, so `summary.txt` says so instead.
"""
function parameterunit(item::_ParamItem)
    key = String(_paramkey(item))
    startswith(key, "log") && (key = key[4:end])
    key == "delta" && return "ppm"
    (startswith(key, "R1_") || startswith(key, "R2_") || key in ("R1", "R2")) &&
        return "s-1"
    (startswith(key, "kex") || startswith(key, "koff") || startswith(key, "kon")) &&
        return "s-1"
    return ""
end

# ---- provenance ---------------------------------------------------------------

"""
    problemcomments(prob) -> Vector{String}

The `#` comment header every CSV carries: the package, the model, the experiments and the
integration region they share.
"""
function problemcomments(prob::ExchangeProblem)
    lines = ["NMRAnalysis.jl $(pkgversion(NMRAnalysis))",
             "Analysis: Exchange (1D)",
             "Model: $(modelname(prob.model))"]
    for expt in prob.experiments
        push!(lines, "Source: $(short_expt_path(expt))")
    end
    push!(lines, "Number of experiments: $(length(prob.experiments))")
    if !isnothing(prob.integration)
        i = prob.integration
        push!(lines,
              "Integration: peak $(i.peakppm) ppm, noise $(i.noiseppm) ppm, " *
              "width $(i.ppmwidth) ppm")
    end
    return lines
end

# ---- the three tables ---------------------------------------------------------

"""
    resultstable(prob) -> (header, rows)

Column names and rows for `results.csv`: one row per experiment, giving its label, its type
and field, how many points it contributed, and the settings constant across it. Keys and no
fitted parameters, every fitted parameter being shared and so global - see the note at the
top of this file.
"""
function resultstable(prob::ExchangeProblem)
    constants = Symbol[]
    for expt in prob.experiments, (name, value) in seriescoordinates(expt)
        value isa AbstractVector || name in constants || push!(constants, name)
    end

    header = ["label", "type", csvcolumn("field", "T"), "points"]
    append!(header, [csvcolumn(name, coordinateunit(name)) for name in constants])

    rows = Vector{String}[]
    for expt in prob.experiments
        coords = Dict(seriescoordinates(expt))
        row = [short_expt_path(expt), experimenttype(expt),
               csvvalue(expt.field_teslas), csvvalue(length(expt.observed_intensities))]
        # A name constant in one experiment may vary in another - nu_SL is a single
        # spin-lock field off-resonance but a whole list of them on-resonance - so where it
        # varies here it reads NA and the values are in series.csv instead.
        for name in constants
            value = get(coords, name, nothing)
            push!(row, value isa AbstractVector ? "NA" : csvvalue(value))
        end
        push!(rows, row)
    end
    return header, rows
end

"""The experiment's own name for itself, as `experimentinfo` reports it."""
function experimenttype(expt::AbstractExperiment)
    for (key, value) in experimentinfo(expt)
        key == "Type" && return value
    end
    return string(nameof(typeof(expt)))
end

"""
    seriestable(prob) -> (header, rows)

Column names and rows for `series.csv`: the measurements, one row per experiment per data
point, with the observed value, its uncertainty and the value the fitted model predicts
there.

Experiments observing different quantities (an intensity, a relaxation rate) get a column
pair each, so no row is ever labelled with a quantity it does not hold.
"""
function seriestable(prob::ExchangeProblem)
    coordnames = Symbol[]
    for expt in prob.experiments, (name, _) in seriescoordinates(expt)
        name in coordnames || push!(coordnames, name)
    end
    observables = unique(observable(expt) for expt in prob.experiments)

    header = ["source", "label"]
    append!(header, [csvcolumn(name, coordinateunit(name)) for name in coordnames])
    for (name, unit) in observables
        append!(header,
                [csvcolumn(name, unit), csvcolumn("$(name)_err", unit),
                 csvcolumn("$(name)_fit", unit)])
    end

    rows = Vector{String}[]
    for expt in prob.experiments
        coords = Dict(seriescoordinates(expt))
        mine = observable(expt)
        label = short_expt_path(expt)
        for i in eachindex(expt.observed_intensities)
            row = [label, label]
            for name in coordnames
                value = get(coords, name, nothing)
                push!(row, csvvalue(value isa AbstractVector ? value[i] : value))
            end
            for obs in observables
                if obs == mine
                    y = expt.observed_intensities[i]
                    append!(row,
                            [csvvalue(Measurements.value(y)),
                             csvvalue(Measurements.uncertainty(y)),
                             csvvalue(expt.predicted_intensities[i])])
                else
                    append!(row, ["NA", "NA", "NA"])
                end
            end
            push!(rows, row)
        end
    end
    return header, rows
end

"""
    globaltable(result) -> (header, rows)

Column names and rows for `global.csv`: every fitted parameter, flat, with the value it
started at, the value it reached, its uncertainty, and whether it was held fixed. This is
where an exchange fit's actual results are, a joint fit having nothing that belongs to one
experiment alone.
"""
function globaltable(result::FitResult)
    header = ["parameter", "value", "error", "unit", "initial", "fixed"]
    initial = Dict(item.flat_index => item
                   for item in _flatten_params_items(result.params0))
    rows = Vector{String}[]
    for item in _flatten_params_items(result.params)
        value = _displayvalue(item, result.params)
        start = haskey(initial, item.flat_index) ?
                _displayvalue(initial[item.flat_index], result.params0) : nothing
        push!(rows,
              [item.label,
               csvvalue(value isa Measurement ? Measurements.value(value) : value),
               csvvalue(value isa Measurement ? Measurements.uncertainty(value) : nothing),
               parameterunit(item),
               csvvalue(start isa Measurement ? Measurements.value(start) : start),
               item.flat_index in result.fixed ? "true" : "false"])
    end
    return header, rows
end

# ---- writing ------------------------------------------------------------------

"""
    writeresults!(result, folder) -> String

Write `results.csv`, `series.csv`, `global.csv` and one `experiments/<label>.csv` per
experiment into `folder`, and return the path of `results.csv`. The per-experiment files
hold that experiment's own rows of `series.csv`, so the data behind each plot sits beside
it under the same basename.
"""
function writeresults!(result::FitResult, folder::AbstractString)
    prob = result.prob
    comments = problemcomments(prob)
    filepath = writetable(joinpath(folder, "results.csv"), comments, resultstable(prob)...)

    header, rows = seriestable(prob)
    writetable(joinpath(folder, "series.csv"), comments, header, rows)

    labelcol = findfirst(==("label"), header)
    for expt in prob.experiments
        label = short_expt_path(expt)
        writetable(joinpath(folder, "experiments", "$(safename(label)).csv"), comments,
                   header, filter(row -> row[labelcol] == label, rows))
    end

    writetable(joinpath(folder, "global.csv"), comments, globaltable(result)...)
    return filepath
end

"""
    writesummary(filepath, result) -> String

Write `summary.txt`: the parameter tables and fit statistics as the REPL shows them, then
the experiment-by-experiment provenance, then the call that repeats the analysis. Replaces
the former pair of `exchange1d_params.txt` and `exchange1d_experiments.txt`, which split
one record across two files.
"""
function writesummary(filepath::AbstractString, result::FitResult)
    backupfile(filepath)
    open(filepath, "w") do io
        show(io, MIME("text/plain"), result)
        println(io)
        writeexperimentsummary(io, result.prob)
        println(io,
                "Concentrations, and any dissociation constant fitted from them, are in " *
                "whatever units the sample metadata gives them in.")
        return nothing
    end
    return filepath
end
