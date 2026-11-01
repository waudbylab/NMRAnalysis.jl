"""
    pre2d(diamagnetic, paramagnetic; type, Trelax=:auto, peaklist=nothing) -> PREExperiment
    pre2d(inputfilenames; concentrations, type, Trelax=:auto, skipplanes=nothing,
          peaklist=nothing) -> PREExperiment

!!! warning "Under development"
    PRE fitting is under active development. Results should be validated carefully before
    use in publication.

Start an interactive GUI for analysing ¹H transverse paramagnetic relaxation enhancement
(PRE) experiments. Every peak is fitted across all the spectra at once, with one
diamagnetic amplitude and linewidths, and the PRE rate Γ₂ accounting for both the line
broadening and the intensity lost during the coherence transfer delays.

Give a diamagnetic and a paramagnetic spectrum for a spin-labelled or metal-tagged protein,
and Γ₂ is in s⁻¹. Give a series of spectra with their `concentrations` for a solvent PRE,
and Γ₂ is in s⁻¹ per unit of concentration, in whatever unit you used.

# Arguments
- `diamagnetic`, `paramagnetic`: processed Bruker data directories, or experiment numbers.
- `inputfilenames`: a list of processed data directories, or experiment numbers, one per
  concentration point.
- `concentrations`: the paramagnetic agent concentration of each spectrum, `0` for the
  diamagnetic reference.
- `type`: `:hsqc` or `:hmqc`. In an HMQC, ¹H is transverse during t₁ as well, so the PRE
  also broadens the indirect dimension.
- `Trelax`: the time (s) for which ¹H magnetisation is transverse outside the acquisition
  and t₁, during which the PRE attenuates the signal. `:auto` takes 1/¹J, the total of the
  transfer delays in a plain HSQC or HMQC, with ¹J = 92 Hz for ¹⁵N and 125 Hz (methyl)
  for ¹³C. Give the time explicitly for a sensitivity-enhanced or other sequence.

# Examples
```julia
# Spin-labelled protein: diamagnetic and paramagnetic samples
pre2d("diamagnetic/pdata/1", "paramagnetic/pdata/1"; type=:hsqc)

# Methyl HMQC, with experiment numbers and an explicit transfer time
pre2d(11, 12; type=:hmqc, Trelax=0.0089)

# Solvent PRE, concentrations in mM
pre2d([10, 11, 12, 13]; concentrations=[0, 1, 5, 10], type=:hmqc)
```

Each peak is plotted as cross-sections through both spectra for a pair, or as its height
against concentration for a series. The summary plots Γ₂, the height in each paramagnetic
spectrum relative to the diamagnetic one, and the diamagnetic linewidths against residue
number. The window blocks until it is closed, and the analysis is returned (see
[`results`](@ref)).
"""
function pre2d(diamagnetic::Union{AbstractString,Integer},
               paramagnetic::Union{AbstractString,Integer}; type, Trelax=:auto,
               peaklist=nothing)
    expt = PREExperiment(asexptpath(diamagnetic), asexptpath(paramagnetic); type, Trelax)
    call = analysiscall("pre2d", diamagnetic, paramagnetic; type, Trelax=expt.Trelax)
    return gui!(expt; peaklist, call)
end

function pre2d(inputfilenames::AbstractVector; concentrations, type, Trelax=:auto,
               skipplanes=nothing, peaklist=nothing)
    skip = checkskipplanes(skipplanes, length(inputfilenames))
    expt = PREExperiment(asexptpath(inputfilenames); concentrations, type, Trelax,
                         skipplanes=skip)
    call = analysiscall("pre2d", inputfilenames; concentrations, type, Trelax=expt.Trelax,
                        skipplanes=nonempty(skip))
    return gui!(expt; peaklist, call)
end

"""
    PREExperiment

Paramagnetic relaxation enhancement experiment.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `concentrations`: Paramagnetic agent concentration of each spectrum, `[0, 1]` for a pair
- `paired`: Whether the spectra are a diamagnetic/paramagnetic pair, so that Γ₂ is in s⁻¹
- `expttype`: Either `:hsqc` or `:hmqc`
- `Trelax`: Transverse ¹H time outside acquisition and t₁ (s)
"""
struct PREExperiment <: ParamagneticExperiment
    specdata::Any
    peaks::Any
    concentrations::Vector{Float64}
    paired::Bool
    expttype::Symbol
    Trelax::Float64

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any
    skipplanes::Vector{Int}

    function PREExperiment(specdata::SpecData, concentrations, paired, expttype, Trelax;
                           skipplanes=Int[])
        checksamples(expttype, concentrations, length(specdata.z), skipplanes)
        Trelax = transfertime(Trelax, first(specdata.nmrdata))
        expt = new(specdata, Observable(Vector{Peak}()), concentrations, paired, expttype,
                   Trelax,
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

primaryparam(::PREExperiment) = :PRE
analysisname(::PREExperiment) = "PRE"

"""
    PREExperiment(diamagnetic, paramagnetic; type, Trelax=:auto)
    PREExperiment(inputfilenames; concentrations, type, Trelax=:auto, skipplanes=Int[])

Create a PRE experiment from a diamagnetic/paramagnetic pair of spectra, or from a series
of spectra at the given paramagnetic agent concentrations. See [`pre2d`](@ref).
"""
function PREExperiment(diamagnetic::AbstractString, paramagnetic::AbstractString; type,
                       Trelax=:auto)
    specdata = preparespecdata([diamagnetic, paramagnetic], ["Diamagnetic", "Paramagnetic"],
                               PREExperiment)
    return PREExperiment(specdata, [0.0, 1.0], true, type, Trelax)
end

function PREExperiment(inputfilenames::AbstractVector; concentrations, type, Trelax=:auto,
                       skipplanes=Int[])
    labels = [c == 0 ? "Diamagnetic" : "Paramagnetic ($c)" for c in concentrations]
    specdata = preparespecdata(inputfilenames, labels, PREExperiment)
    return PREExperiment(specdata, collect(Float64, concentrations), false, type, Trelax;
                         skipplanes)
end

# ---- the model -------------------------------------------------------------------

function modelparameters(expt::PREExperiment)
    cmax = maximum(expt.concentrations)
    return [(key=:PRE, label="Γ₂", initial=10.0 / cmax, lower=0.0, upper=1000.0 / cmax,
             scale=1 / cmax)]
end

planecondition(expt::PREExperiment, i) = (c=expt.concentrations[i],)

function peakmodel(expt::PREExperiment, cond, R2x, R2y, extras)
    return transverseplane(expt, cond.c, R2x, R2y, extras[1])
end

function referenceplanes(expt::PREExperiment, i)
    return filter(j -> expt.concentrations[j] == 0,
                  usedplanes(expt))
end
function planelabel(expt::PREExperiment, i)
    return expt.paired ? "Paramagnetic" :
           "c = $(expt.concentrations[i])"
end

function summarypanels(expt::PREExperiment, peaks)
    return [ratepanel(expt, peaks, :PRE, "Γ₂"), intensitypanel(expt, peaks),
            transversepanel(expt, peaks)]
end

# ---- the peak plot of a series: height against concentration -----------------------

curvelabels(::PREExperiment) = (xlabel="Concentration", ylabel="Relative peak height")

function drawcurve!(ax, peak, expt::PREExperiment)
    c = expt.concentrations
    A = peak.parameters[:amp0].value[][1]
    y = peak.parameters[:amp].value[] ./ A
    e = peak.parameters[:amp].uncertainty[] ./ A
    used = [i ∉ skipset(expt) for i in eachindex(c)]
    errorbars!(ax, c[used], y[used], e[used]; whiskerwidth=6, color=seriescolour(1))
    scatter!(ax, c[used], y[used]; color=seriescolour(1))
    any(!, used) &&
        scatter!(ax, c[.!used], y[.!used]; color=:transparent, strokecolor=:gray60,
                 strokewidth=1.5)
    peak.postfitted[] || return nothing
    φp = fitvalues(peak, expt)
    cs = range(0, 1.05 * maximum(c), 100)
    return lines!(ax, cs,
                  [peakmodel(expt, (c=ci,), φp[3], φp[4], φp[5:end])[3] for ci in cs];
                  color=:red)
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::PREExperiment)
    samples = expt.paired ? "diamagnetic, paramagnetic" :
              "concentrations $(join(expt.concentrations, ", "))"
    return "Analysis type: PRE experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Samples: $samples\n" *
           "Experiment type: $(expt.expttype == :hsqc ? "HSQC" : "HMQC")\n" *
           "Transfer time Trelax / ms: $(round(1000 * expt.Trelax; digits=2))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end
