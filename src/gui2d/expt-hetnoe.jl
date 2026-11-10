"""
    hetnoe2d(reference, saturated; skipplanes=nothing, peaklist=nothing) -> HetNOEExperiment
    hetnoe2d(planes; saturation=nothing, skipplanes=nothing, peaklist=nothing,
             prompt=isinteractive()) -> HetNOEExperiment

Start interactive GUI for analysing 2D heteronuclear NOE data.

Give the `reference` and `saturated` spectra, each a single spectrum or a list of them, or
give every plane in `planes` with `saturation` saying which are saturated. Spectra can be
paths or Bruker experiment numbers, individually or as a list or range. With several
pairs, results are averaged across them.

`saturation` is one `true` or `false` per plane. Without it, it is read from each
spectrum's `hetnoe.saturated` annotation, and otherwise you are asked.

`skipplanes` lists planes (1-based, in the order of `planes`, or with reference and
saturated planes interleaved) left out of the averages; at least one of each must remain.
`peaklist` is a peak list to start from, or a folder holding one. The window blocks until it
is closed, and the analysis is returned (see [`results`](@ref)).

# Examples
```julia
hetnoe2d("expno1/pdata/231", "expno1/pdata/232")      # one reference / saturated pair
hetnoe2d(["1/pdata/1", "2/pdata/1"], ["3/pdata/1", "4/pdata/1"])   # two pairs, averaged
hetnoe2d(1, 2)                                         # Bruker experiment numbers
hetnoe2d(10:13; saturation=[false, true, false, true]) # planes, in any order
```
"""
function hetnoe2d(reference, saturated; kwargs...)
    ref, sat = planelist(reference), planelist(saturated)
    length(ref) == length(sat) ||
        throw(ArgumentError("got $(length(ref)) reference and $(length(sat)) saturated " *
                            "spectra: give one of each per pair"))
    planes = collect(Iterators.flatten(zip(ref, sat)))
    return hetnoe2d(planes; saturation=repeat([false, true], length(ref)), kwargs...)
end

function hetnoe2d(planes; saturation=nothing, skipplanes=nothing, peaklist=nothing,
                  prompt::Bool=isinteractive())
    planes = planelist(planes)
    saturation = @something(saturation, annotatedsaturation(planes),
                            asksaturation(length(planes); prompt))
    length(saturation) == length(planes) ||
        throw(ArgumentError("got $(length(saturation)) saturation values for " *
                            "$(length(planes)) planes"))
    skip = checkskipplanes(skipplanes, length(saturation))
    expt = HetNOEExperiment(planes, collect(Bool, saturation); skipplanes=skip)
    call = analysiscall("hetnoe2d", planes; saturation=collect(Bool, saturation),
                        skipplanes=nonempty(skip))
    return gui!(expt; peaklist, call)
end

"""Which planes are saturated, from the `hetnoe.saturated` annotations, or `nothing`."""
function annotatedsaturation(planes)
    values = [annotation(loadnmr(p), :hetnoe, :saturated) for p in planes]
    any(isnothing, values) && return nothing
    flat = reduce(vcat, values)
    return length(flat) == length(planes) ? Bool.(flat) : nothing
end

function asksaturation(n; prompt::Bool=true)
    prompt || throw(ArgumentError("which planes are saturated could not be read from the " *
                                  "annotations: pass saturation, e.g. [false, true]"))
    values = askvector("saturation flags (1 for saturated, 0 for reference)", n)
    return values .!= 0
end

"""
    HetNOEExperiment

Heteronuclear NOE experiment with reference and saturated spectra.

# Fields
- `specdata`: Spectral data and metadata  
- `peaks`: Observable list of peaks
- `saturation`: Vector of booleans indicating saturated spectra
"""
struct HetNOEExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    saturation::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any
    skipplanes::Vector{Int}

    function HetNOEExperiment(specdata, peaks, saturation; skipplanes=Int[])
        used = [i ∉ skipplanes for i in eachindex(saturation)]
        any(saturation .& used) && any(.!saturation .& used) ||
            throw(ArgumentError("skipplanes must leave a reference and a saturated plane"))
        expt = new(specdata, peaks, saturation,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}(),
                   collect(Int, skipplanes))
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

struct HetNOEVisualisation <: VisualisationStrategy end
visualisationtype(::HetNOEExperiment) = HetNOEVisualisation()
primaryparam(::HetNOEExperiment) = :hetnoe

"""
    HetNOEExperiment(planefilenames, saturation)

Create hetNOE experiment from a list of input planes and a list of
true/false values indicating where saturation has been applied.
"""
function HetNOEExperiment(planefilenames, saturation::AbstractVector{Bool};
                          skipplanes=Int[])
    specdata = preparespecdata(planefilenames, saturation, HetNOEExperiment)
    peaks = Observable(Vector{Peak}())

    return HetNOEExperiment(specdata, peaks, collect(saturation); skipplanes)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(planefilenames, saturation, ::Type{HetNOEExperiment})
    @debug "Preparing spec data for hetNOE experiment: $planefilenames"
    spectra = map(loadnmr, planefilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(sat -> sat ? "Sat" : "Ref", saturation)

    return SpecData(spectra, x, y,
                    z ./ σ[1],
                    σ ./ σ[1],
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::HetNOEExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.0)
    R2y0 = MaybeVector(10.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    newpeak.postparameters[:hetnoe] = Parameter("Heteronuclear NOE", 0.7)
    newpeak.postparameters[:amp] = Parameter("Amplitude", maximum(amp0))

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::HetNOEExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog = 10
    sat = expt.saturation
    used = [i ∉ skipset(expt) for i in eachindex(sat)]
    A = peak.parameters[:amp].value[] .± peak.parameters[:amp].uncertainty[]

    Iref = A[.!sat .& used]
    Isat = A[sat .& used]
    Iref = sum(Iref) / length(Iref)
    Isat = sum(Isat) / length(Isat)
    hetNOE = Isat / Iref
    peak.postparameters[:hetnoe].value[] .= Measurements.value(hetNOE)
    peak.postparameters[:hetnoe].uncertainty[] .= Measurements.uncertainty(hetNOE)
    peak.postparameters[:amp].value[] .= Measurements.value(Iref)
    peak.postparameters[:amp].uncertainty[] .= Measurements.uncertainty(Iref)

    return peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::HetNOEExperiment, idx)
    if expt.saturation[idx]
        "Saturated ($idx of $(nslices(expt)))"
    else
        "Reference ($idx of $(nslices(expt)))"
    end
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::HetNOEExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "HetNOE: $(peak.postparameters[:hetnoe].value[][1] ± peak.postparameters[:hetnoe].uncertainty[][1])\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "Amplitude: $(peak.postparameters[:amp].value[][1] ± peak.postparameters[:amp].uncertainty[][1])\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
end

"""Return formatted text describing experiment."""
analysisname(::HetNOEExperiment) = "heteronuclear NOE"

function experimentinfo(expt::HetNOEExperiment)
    return "Analysis type: Heteronuclear NOE\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Saturation: $(join(expt.saturation, ", "))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

## visualisation
function get_hetnoe_data(peak, expt::HetNOEExperiment)
    x = 1:nslices(expt)

    if isnothing(peak)
        y = fill(0.0, nslices(expt))
        err = [(1.0 * i, 0.0, 0.0) for i in 1:nslices(expt)]
        fit = [0.0]
        return (x, y, err, fit)
    end

    Iref = peak.postparameters[:amp].value[][1]
    hetnoe = peak.postparameters[:hetnoe].value[][1]

    y = [peak.parameters[:amp].value[][i] / Iref for i in 1:nslices(expt)]
    err = [(1.0 * i,
            peak.parameters[:amp].value[][i] / Iref,
            peak.parameters[:amp].uncertainty[][i] / Iref)
           for i in 1:nslices(expt)]
    fit = [1, hetnoe]

    return (x, y, err, fit)
end

function completestate!(state, expt, ::HetNOEVisualisation)
    @debug "completing state for hetNOE visualisation"
    state[:peak_plot_data] = lift(peak -> get_hetnoe_data(peak, expt), state[:current_peak])
    state[:peak_plot_x] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_y] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_err] = lift(d -> d[3], state[:peak_plot_data])
    return state[:peak_plot_fit] = lift(d -> d[4], state[:peak_plot_data])
end

function plot_peak!(panel, peak, expt, ::HetNOEVisualisation)
    x, y, err, fit = get_hetnoe_data(peak, expt)

    ax = Axis(panel[1, 1];
              xlabel="Experiment",
              ylabel="Relative amplitude",
              xticks=(1:nslices(expt), expt.specdata.zlabels))

    hlines!(ax, [0]; color=:grey50, linewidth=1)
    hlines!(ax, fit; linewidth=2, color=:red)
    errorbars!(ax, err; whiskerwidth=10)
    return barplot!(ax, x, y)
end

function makepeakplot!(gui, state, expt, ::HetNOEVisualisation)
    @debug "making peak plot for hetNOE visualisation"
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel="Experiment",
                                 ylabel="Relative amplitude",
                                 xticks=(1:nslices(expt), expt.specdata.zlabels))

    hlines!(ax, [0]; color=:grey50, linewidth=1)
    hlines!(ax, state[:peak_plot_fit]; linewidth=2, color=:red)
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    return barplot!(ax, state[:peak_plot_x], state[:peak_plot_y])
end
