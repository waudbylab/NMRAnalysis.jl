# Output tables, written to the layout and column rules in
# `docs/src/advanced/conventions.md`:
#
#   out/
#     summary.txt        the record to read
#     results.csv        one row per peak: positions, linewidths and derived parameters
#     series.csv         the measurements, one row per peak per plane
#     global.csv         parameters fitted once across every peak (a titration Kd)
#     summary.pdf
#     peaks/<label>.pdf, peaks/<label>.csv
#     cluster_*.pdf
#
# The **entity** is the peak. Anything that varies plane by plane - the amplitude always,
# and for a moving-peak experiment the positions and linewidths too - belongs in
# `series.csv`, keyed by the plane's own coordinate rather than by an `amp[7]` column whose
# meaning lived only in a comment line.

# ---- coordinates --------------------------------------------------------------

"""
    seriescoordinates(expt) -> Vector{Pair{Symbol,Any}}

What distinguishes one plane from another, as `name => values`. A vector gives one entry
per plane; a scalar is a setting constant across the experiment (a saturation time, a CPMG
relaxation delay) and is repeated down the rows.

Named per quantity rather than reduced to one `x` column, so that a saved file says what
its planes actually were. This is the same hook Exchange1D defines for the same reason.

The two generic experiment types differ only by the model fitted, so they delegate to a
second method dispatching on that: `rdc2d` needs two keys where a titration needs one.
"""
seriescoordinates(e::IntensityExperiment) = seriescoordinates(e, e.model)
seriescoordinates(e::MovingExperiment) = seriescoordinates(e, e.model)
seriescoordinates(e, model::FittingModel) = [coordinatename(model) => e.x]

# Nothing was arrayed, so the plane index is the coordinate - and `series.csv` already
# carries it as a column of its own.
seriescoordinates(e, ::NoFitting) = Pair{Symbol,Any}[]

function seriescoordinates(e::CESTExperiment)
    return [:offset => e.frequencies, :B1 => e.B1, :Tsat => e.Tsat]
end
seriescoordinates(e::CPMGExperiment) = [:nu_CPMG => e.vCPMG, :Trelax => e.Trelax]
seriescoordinates(e::HetNOEExperiment) = [:saturated => collect(e.saturation)]
seriescoordinates(e::CCRExperiment) = [:buildup => collect(e.isbuildup), :Trelax => e.T]

"""
    coordinatename(model) -> Symbol

What the per-plane `x` of an experiment fitted with `model` actually is. The generic
fallback is the uninformative `:x`, which is the honest answer for a model fitted against
whatever the caller supplied (`modelfit2d`).
"""
coordinatename(::FittingModel) = :x
coordinatename(::ExponentialModel) = :time
coordinatename(::RecoveryModel) = :time
coordinatename(::MethylCCRModel) = :time
coordinatename(::TitrationModel) = :concentration

const COORDINATE_UNITS = Dict(:time => "s", :offset => "ppm", :B1 => "Hz", :Tsat => "s",
                              :nu_CPMG => "Hz", :Trelax => "s")

coordinateunit(name::Symbol) = get(COORDINATE_UNITS, name, "")

# ---- units --------------------------------------------------------------------

# ASCII units for the parameters shared by more than one experiment, as they appear in
# column headers. A parameter missing here gets no unit rather than a guessed one; `Kd` is
# the notable blank, a titration's concentrations coming from the sample metadata, so its
# units are whatever that metadata used.
const PARAM_UNITS = Dict(:x => "ppm", :y => "ppm", :R2x => "s-1", :R2y => "s-1",
                         :R => "s-1", :R1 => "s-1", :R2 => "s-1", :R20 => "s-1",
                         :PRE => "s-1", :eta => "s-1", :S2tc => "ns",
                         :CSP => "ppm", :dX => "ppm", :dY => "ppm",
                         :Xfree => "ppm", :Xbound => "ppm",
                         :Yfree => "ppm", :Ybound => "ppm")

"""
    paramunit(expt, name) -> String

ASCII unit for a fitted or derived parameter, as it appears in a column header.
Experiment-dispatched, so an experiment introducing its own parameters declares their
units in its own `expt-*.jl` rather than here; the default falls back to the shared table
above and then to no unit.
"""
paramunit(::Experiment, name::Symbol) = get(PARAM_UNITS, name, "")
paramunit(name::Symbol) = get(PARAM_UNITS, name, "")

"""
    globalparams(expt) -> Vector{Symbol}

Derived parameters fitted once across every peak rather than per peak, and so belonging in
`global.csv`. A titration's `Kd` is fitted globally and then copied onto every peak, so
writing it in `results.csv` would repeat one number down a column.

Declared beside the `postfitglobal!` that computes it - see `globalparams(::TitrationModel)`
in `expt-moving.jl`.
"""
globalparams(::Experiment) = Symbol[]
globalparams(::FittingModel) = Symbol[]
globalparams(e::IntensityExperiment) = globalparams(e.model)
globalparams(e::MovingExperiment) = globalparams(e.model)

"""
    fittedamplitudes(peak, expt) -> Vector{Float64}

The model evaluated at each plane's own coordinate, for the `amp_fit` column of
`series.csv`, so that a residual is a subtraction. `NaN` (written `NA`) for every plane
where nothing was fitted through the amplitudes themselves: a heteronuclear NOE, a CCR
rate and a CEST profile are fitted from ratios or from transformed intensities, so no
curve passes through the amplitudes to report.
"""
fittedamplitudes(peak, expt) = fill(NaN, nslices(expt))

# ---- provenance ---------------------------------------------------------------

"""Comment lines heading every CSV: the experiment description plus the fitting radii,
which `readpeaklist!` reads back."""
function resultcomments(expt)
    lines = ["NMRAnalysis.jl $(pkgversion(GUI2D))"]
    append!(lines, split(experimentinfo(expt), '\n'))
    push!(lines, "X radius / ppm: $(round(expt.xradius[]; digits=4))")
    push!(lines, "Y radius / ppm: $(round(expt.yradius[]; digits=4))")
    return lines
end

"""Where plane `i` came from. One file per plane for a multi-file experiment; the same
pseudo-3D dataset for every plane otherwise."""
function planesource(expt, i)
    nmrdata = expt.specdata.nmrdata
    spec = nmrdata[min(i, length(nmrdata))]
    return string(something(spec[:filename], ""))
end

# ---- results.csv --------------------------------------------------------------

"The lineshape parameters of a peak: where it is, and how broad it is in each dimension."
const POSITION_PARAMS = (:x, :y, :R2x, :R2y)

"""
    derivedkeys(expt) -> Vector{Symbol}

Derived parameter names for `results.csv`, the experiment's [`primaryparam`](@ref) first
and anything global left out.
"""
function derivedkeys(expt)
    peaks = expt.peaks[]
    isempty(peaks) && return Symbol[]
    keys_ = [k for k in keys(first(peaks).postparameters) if !(k in globalparams(expt))]
    primary = primaryparam(expt)
    i = findfirst(==(primary), keys_)
    isnothing(i) || pushfirst!(keys_, popat!(keys_, i))
    return keys_
end

"""
    resultstable(expt) -> (header, rows)

Column names and rows for `results.csv`: one row per peak, its identity and the parameters
derived from the fit. This is the table to plot against residue number.

Positions and linewidths are not here. They are per-plane quantities that happen to take
one value per plane when the peaks are fixed, so they live in `series.csv` whether or not
they move, and the layout is the same either way.
"""
function resultstable(expt)
    derived = derivedkeys(expt)

    header = ["label", "resnum", "resname", "atom"]
    for k in derived
        append!(header, collect(csvcolumns(k, paramunit(expt, k))))
    end

    rows = Vector{String}[]
    for peak in sortedpeaks(expt)
        lbl = parse_label(peak.label[])
        row = [peak.label[], string(lbl.resnum),
               lbl.onelettercode == '?' ? "" : string(lbl.onelettercode), lbl.atom]
        for k in derived
            push!(row, format_post(peak, k, :value))
            push!(row, format_post(peak, k, :uncertainty))
        end
        push!(rows, row)
    end
    return header, rows
end

# ---- peaklist.csv: the user's input -------------------------------------------

"""
    peaklisttable(expt) -> (header, rows)

Column names and rows for `peaklist.csv`: where the user *put* each peak, as opposed to
where the fit moved it to. One row per peak, or one row per peak per plane where the
positions were placed plane by plane.

The distinction matters because these are different things with different lives. A peak
list is an input: you pick it once, reuse it on another dataset, hand it to a colleague, or
import it from elsewhere. A fitted position is an output of one particular fit. Reading
back the fitted positions as the next run's starting point - which is what loading
`results.csv` used to do - conflates the two, and for a moving-peak experiment there is no
single fitted position to read.

`plane` is left blank for a peak whose initial position is one value for the whole
experiment, meaning the row applies to every plane, and carries the plane number where the
peak was tracked across them. The radii are per-peak and repeat down that peak's rows.
"""
function peaklisttable(expt)
    header = ["label", "plane", csvcolumn("x", "ppm"), csvcolumn("y", "ppm"),
              csvcolumn("xradius", "ppm"), csvcolumn("yradius", "ppm")]
    rows = Vector{String}[]
    for peak in sortedpeaks(expt)
        x0 = peak.parameters[:x].initialvalue[]
        y0 = peak.parameters[:y].initialvalue[]
        rx, ry = csvvalue(peak.xradius[]), csvvalue(peak.yradius[])
        # A `SingleElementVector` holds one position for the whole experiment; a peak
        # tracked plane by plane holds one per plane.
        if x0 isa SingleElementVector && y0 isa SingleElementVector
            push!(rows, [peak.label[], "", csvvalue(x0[1]), csvvalue(y0[1]), rx, ry])
        else
            for i in 1:nslices(expt)
                push!(rows,
                      [peak.label[], string(i), csvvalue(x0[i]), csvvalue(y0[i]), rx, ry])
            end
        end
    end
    return header, rows
end

# ---- series.csv ---------------------------------------------------------------

"""
    seriestable(expt) -> (header, rows)

Column names and rows for `series.csv`: one row per peak per plane, with the plane's
coordinates and everything measured there - the amplitude, the position and the linewidths.
A fixed-peak experiment simply repeats its one position down the rows.

`plane` is always written and is load-bearing where the planes share one file: every row of
a pseudo-3D experiment has the same `source`, and the plane index is then the only thing
telling two rows apart.
"""
function seriestable(expt)
    n = nslices(expt)
    coords = seriescoordinates(expt)

    header = ["source", "label", "plane"]
    append!(header, [csvcolumn(name, coordinateunit(name)) for (name, _) in coords])
    append!(header, collect(csvcolumns(:amp, paramunit(expt, :amp))))
    push!(header, csvcolumn("amp_fit", paramunit(expt, :amp)))
    for k in POSITION_PARAMS
        append!(header, collect(csvcolumns(k, paramunit(expt, k))))
    end

    rows = Vector{String}[]
    for peak in sortedpeaks(expt)
        ampfit = fittedamplitudes(peak, expt)
        for i in 1:n
            row = [planesource(expt, i), peak.label[], string(i)]
            for (_, value) in coords
                push!(row, csvvalue(value isa AbstractVector ? value[i] : value))
            end
            push!(row, format_param(peak, :amp, i, :value))
            push!(row, format_param(peak, :amp, i, :uncertainty))
            push!(row, csvvalue(ampfit[i]))
            for k in POSITION_PARAMS
                push!(row, format_param(peak, k, i, :value))
                push!(row, format_param(peak, k, i, :uncertainty))
            end
            push!(rows, row)
        end
    end
    return header, rows
end

# ---- global.csv ---------------------------------------------------------------

"""
    globaltable(expt) -> (header, rows)

Column names and rows for `global.csv`: the parameters [`globalparams`](@ref) names, read
from the first peak that carries them, since a global fit gives every peak the same value.
"""
function globaltable(expt)
    header = ["parameter", "value", "error", "unit"]
    rows = Vector{String}[]
    peaks = expt.peaks[]
    for k in globalparams(expt)
        i = findfirst(p -> haskey(p.postparameters, k), peaks)
        isnothing(i) && continue
        push!(rows,
              [string(k), format_post(peaks[i], k, :value),
               format_post(peaks[i], k, :uncertainty), paramunit(expt, k)])
    end
    return header, rows
end

# ---- writing ------------------------------------------------------------------

"""
    writeresults!(expt, folder) -> String

Write `peaklist.csv` (what was picked), `series.csv`, `results.csv` and `global.csv` (each
where there is anything to put in it) and one `peaks/<label>.csv` per peak into `folder`,
and return the path of the series file. The per-peak files hold that peak's own rows of
`series.csv`, so the data behind each plot sits beside it under the same basename.

`results.csv` is skipped by an experiment that derives nothing per peak (`fit2d`,
`peaktrack2d`): it would carry peak labels and no values.
"""
function writeresults!(expt, folder)
    comments = resultcomments(expt)

    # What the user picked, kept apart from what the fit produced - see `peaklisttable`.
    writetable(joinpath(folder, "peaklist.csv"), comments, peaklisttable(expt)...)

    isempty(derivedkeys(expt)) ||
        writetable(joinpath(folder, "results.csv"), comments, resultstable(expt)...)

    header, rows = seriestable(expt)
    filepath = writetable(joinpath(folder, "series.csv"), comments, header, rows)

    labelcol = findfirst(==("label"), header)
    for peak in expt.peaks[]
        label = peak.label[]
        writetable(joinpath(folder, "peaks", "$(safename(label)).csv"), comments, header,
                   filter(row -> row[labelcol] == label, rows))
    end

    gheader, grows = globaltable(expt)
    isempty(grows) || writetable(joinpath(folder, "global.csv"), comments, gheader, grows)
    return filepath
end

"""
    writesummary(filepath, expt) -> String

Write `summary.txt`: the experiment description, the fitting radii, and the headline
parameter for every peak, rounded for reading rather than written at full precision.
"""
function writesummary(filepath, expt)
    backupfile(filepath)
    peaks = sortedpeaks(expt)
    primary = primaryparam(expt)
    open(filepath, "w") do f
        for line in resultcomments(expt)
            isempty(strip(line)) && continue
            println(f, line)
        end
        println(f)
        gheader, grows = globaltable(expt)
        if !isempty(grows)
            println(f, "Global parameters:")
            for row in grows
                println(f, "  $(row[1]): $(row[2]) +/- $(row[3]) $(row[4])")
            end
            println(f)
        end
        # An experiment that derives nothing per peak (fit2d, peaktrack2d) has :amp as its
        # primary parameter, which is per plane and so lives in `parameters`. Say where the
        # numbers are rather than print a heading with nothing under it.
        if !any(haskey(p.postparameters, primary) for p in peaks)
            println(f, "Nothing is derived per peak. The fitted amplitudes, positions and")
            println(f, "linewidths of each peak in each plane are in series.csv.")
            return nothing
        end
        unit = paramunit(expt, primary)
        println(f, "$(primary)$(isempty(unit) ? "" : " / $unit") by peak:")
        for peak in peaks
            haskey(peak.postparameters, primary) || continue
            value = tryparse(Float64, format_post(peak, primary, :value))
            err = tryparse(Float64, format_post(peak, primary, :uncertainty))
            isnothing(value) && continue
            rounded = round(value; sigdigits=4)
            suffix = isnothing(err) ? "" : " +/- $(round(err; sigdigits=2))"
            println(f, "  $(rpad(peak.label[], 12)) $rounded$suffix")
        end
        return nothing
    end
    return filepath
end
