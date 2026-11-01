"""
    lpre2d(diamagnetic, paramagnetic; type, recycledelays=nothing, Trelax=:auto,
           peaklist=nothing, prompt=isinteractive()) -> LPREExperiment
    lpre2d(inputfilenames; concentrations, type, recycledelays=nothing, Trelax=:auto,
           skipplanes=nothing, peaklist=nothing, prompt=isinteractive()) -> LPREExperiment

!!! warning "Under development"
    PRE fitting is under active development. Results should be validated carefully before
    use in publication.

Start an interactive GUI for analysing longitudinal and transverse ¹H paramagnetic
relaxation enhancements together, from diamagnetic and paramagnetic spectra recorded at a
series of recycle delays. Every peak is fitted across all the spectra at once, with one
diamagnetic amplitude, linewidths and R₁, and the enhancements Γ₁ and Γ₂. Γ₂ is then free of
the bias that incomplete recovery gives at short delays.

Give lists of diamagnetic and paramagnetic spectra for a spin-labelled or metal-tagged
protein, and Γ₁ and Γ₂ are in s⁻¹. Give one list with its `concentrations` for a solvent
PRE, and they are in s⁻¹ per unit of concentration, in whatever unit you used.

# Arguments
- `diamagnetic`, `paramagnetic`: lists of processed Bruker data directories, or experiment
  numbers, one per recycle delay.
- `inputfilenames`: a list of processed data directories, or experiment numbers.
- `concentrations`: the paramagnetic agent concentration of each spectrum, `0` for the
  diamagnetic reference.
- `type`: `:hsqc` or `:hmqc`, as for [`pre2d`](@ref).
- `recycledelays`: the time (s) the ¹H magnetisation recovers between scans, one per
  spectrum (diamagnetic first), or for two lists one list applied to both. Used exactly as
  given. By default each spectrum's d1 plus its acquisition time.
- `Trelax`: the transfer time, as for [`pre2d`](@ref).

# Examples
```julia
# Spin-labelled protein, nine recycle delays each
lpre2d(10:2:26, 11:2:27; type=:hmqc)

# Solvent PRE with 2 mM Gd, nine delays without and nine with
lpre2d([10:18; 103:111]; concentrations=[fill(0, 9); fill(2, 9)], type=:hmqc)
```

Each peak is plotted as its height against recovery time, a series and fitted curve for
each concentration. The summary plots Γ₂, Γ₁, the height in each paramagnetic spectrum
relative to the diamagnetic one at the same delay, and the diamagnetic R₂ and R₁ against
residue number. The window blocks until it is closed, and the analysis is returned (see
[`results`](@ref)).
"""
function lpre2d(diamagnetic::AbstractVector, paramagnetic::AbstractVector; type,
                recycledelays=nothing, Trelax=:auto, peaklist=nothing,
                prompt::Bool=isinteractive())
    expt = LPREExperiment(asexptpath(diamagnetic), asexptpath(paramagnetic); type,
                          recycledelays, Trelax, prompt)
    call = analysiscall("lpre2d", diamagnetic, paramagnetic; type,
                        recycledelays=expt.recoverytimes, Trelax=expt.Trelax)
    return gui!(expt; peaklist, call)
end

function lpre2d(inputfilenames::AbstractVector; concentrations, type, recycledelays=nothing,
                Trelax=:auto, skipplanes=nothing, peaklist=nothing,
                prompt::Bool=isinteractive())
    skip = checkskipplanes(skipplanes, length(inputfilenames))
    expt = LPREExperiment(asexptpath(inputfilenames); concentrations, type, recycledelays,
                          Trelax, skipplanes=skip, prompt)
    call = analysiscall("lpre2d", inputfilenames; concentrations, type,
                        recycledelays=expt.recoverytimes, Trelax=expt.Trelax,
                        skipplanes=nonempty(skip))
    return gui!(expt; peaklist, call)
end

"""
    LPREExperiment

Longitudinal and transverse paramagnetic relaxation enhancement experiment.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `concentrations`: Paramagnetic agent concentration of each spectrum, 0 or 1 for a pair
- `paired`: Whether the spectra are diamagnetic/paramagnetic pairs, so that Γ is in s⁻¹
- `expttype`: Either `:hsqc` or `:hmqc`
- `Trelax`: Transverse ¹H time outside acquisition and t₁ (s)
- `recoverytimes`: Time the ¹H magnetisation recovers before each scan (s)
"""
struct LPREExperiment <: ParamagneticExperiment
    specdata::Any
    peaks::Any
    concentrations::Vector{Float64}
    paired::Bool
    expttype::Symbol
    Trelax::Float64
    recoverytimes::Vector{Float64}

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any
    skipplanes::Vector{Int}

    function LPREExperiment(specdata::SpecData, concentrations, paired, expttype, Trelax,
                            recoverytimes; skipplanes=Int[])
        n = length(specdata.z)
        checksamples(expttype, concentrations, n, skipplanes)
        length(recoverytimes) == n ||
            throw(ArgumentError("got $(length(recoverytimes)) recovery times for $n spectra"))
        all(>(0), recoverytimes) ||
            throw(ArgumentError("recovery times must be positive"))
        length(unique(recoverytimes[i] for i in 1:n if i ∉ skipplanes)) ≥ 2 ||
            throw(ArgumentError("LPRE needs spectra at two or more recycle delays"))
        Trelax = transfertime(Trelax, first(specdata.nmrdata))
        expt = new(specdata, Observable(Vector{Peak}()), concentrations, paired, expttype,
                   Trelax, recoverytimes,
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

primaryparam(::LPREExperiment) = :Gamma2
analysisname(::LPREExperiment) = "longitudinal PRE"

"""
    LPREExperiment(diamagnetic, paramagnetic; type, recycledelays=nothing, Trelax=:auto,
                   prompt=isinteractive())
    LPREExperiment(inputfilenames; concentrations, type, recycledelays=nothing,
                   Trelax=:auto, skipplanes=Int[], prompt=isinteractive())

Create an LPRE experiment from lists of diamagnetic and paramagnetic spectra, or from one
list at the given paramagnetic agent concentrations. See [`lpre2d`](@ref).
"""
function LPREExperiment(diamagnetic::AbstractVector, paramagnetic::AbstractVector; type,
                        recycledelays=nothing, Trelax=:auto, prompt=isinteractive())
    nd, np = length(diamagnetic), length(paramagnetic)
    if !isnothing(recycledelays) && length(recycledelays) == nd == np
        recycledelays = [recycledelays; recycledelays]
    end
    return loadlpre([diamagnetic; paramagnetic], [zeros(nd); ones(np)], true, type, Trelax,
                    recycledelays, Int[], prompt)
end

function LPREExperiment(inputfilenames::AbstractVector; concentrations, type,
                        recycledelays=nothing, Trelax=:auto, skipplanes=Int[],
                        prompt=isinteractive())
    return loadlpre(inputfilenames, collect(Float64, concentrations), false, type, Trelax,
                    recycledelays, skipplanes, prompt)
end

function loadlpre(files, concentrations, paired, type, Trelax, recycledelays, skipplanes,
                  prompt)
    checksamples(type, concentrations, length(files), skipplanes)
    specdata = preparespecdata(files, fill("", length(files)), LPREExperiment)
    τ = recoverytimes(specdata.nmrdata, recycledelays; prompt)
    for (i, c) in enumerate(concentrations)
        sample = c == 0 ? "Diamagnetic" : paired ? "Paramagnetic" : "Paramagnetic ($c)"
        specdata.zlabels[i] = "$sample, τ = $(round(τ[i]; sigdigits=3)) s"
    end
    return LPREExperiment(specdata, concentrations, paired, type, Trelax, τ; skipplanes)
end

"""
    recoverytimes(spectra, recycledelays; prompt) -> Vector{Float64}

`recycledelays` exactly as given, or else each spectrum's d1 plus its acquisition time,
TD/(2 SW_h) in the directly detected dimension, during which the ¹H magnetisation also
recovers. Asked for where the acquisition parameters are missing.
"""
function recoverytimes(spectra, recycledelays; prompt)
    n = length(spectra)
    if !isnothing(recycledelays)
        τ = collect(Float64, recycledelays)
        length(τ) == n ||
            throw(ArgumentError("got $(length(τ)) recycle delays for $n spectra"))
        return τ
    end
    τ = map(spectra) do spec
        d1, td, swh = acqusvalue(spec, :d, 1), acqusvalue(spec, :td), acqusvalue(spec, :sw_h)
        return any(isnothing, (d1, td, swh)) ? nothing : d1 + td / (2 * swh)
    end
    any(isnothing, τ) || return Vector{Float64}(τ)
    return askvector("recovery times (d1 plus acquisition time)", n; unit="s", prompt)
end

# ---- the model -------------------------------------------------------------------

function modelparameters(expt::LPREExperiment)
    cmax = maximum(expt.concentrations)
    return [(key=:R1, label="R₁", initial=1.5, lower=0.05, upper=50.0, scale=1.0),
            (key=:Gamma1, label="Γ₁", initial=1.0 / cmax, lower=0.0, upper=100.0 / cmax,
             scale=1 / cmax),
            (key=:Gamma2, label="Γ₂", initial=10.0 / cmax, lower=0.0,
             upper=1000.0 / cmax, scale=1 / cmax)]
end

planecondition(expt::LPREExperiment, i) = (c=expt.concentrations[i],
                                           τ=expt.recoverytimes[i])

"""
    peakmodel(expt::LPREExperiment, cond, R2x, R2y, (R1, Γ1, Γ2)) -> (R2x, R2y, height)

The transverse model of `transverseplane`, with the height scaled by the recovery of the
¹H magnetisation from saturation, 1 - exp(-(R₁ + Γ₁c)τ), after a recovery time τ
(Pintacuda & Otting, J. Am. Chem. Soc. 124, 372 (2002)). R₁ and Γ₁ are therefore apparent
rates for the sequence used, which a selective or band-selective sequence makes faster.
"""
function peakmodel(expt::LPREExperiment, cond, R2x, R2y, extras)
    R1, Γ1, Γ2 = extras
    R2xp, R2yp, s = transverseplane(expt, cond.c, R2x, R2y, Γ2)
    return R2xp, R2yp, s * (1 - exp(-(R1 + Γ1 * cond.c) * cond.τ))
end

# a paramagnetic plane is compared with the diamagnetic planes at the same delay
function referenceplanes(expt::LPREExperiment, i)
    τ = expt.recoverytimes
    return filter(j -> expt.concentrations[j] == 0 && isapprox(τ[j], τ[i]; rtol=1e-6),
                  usedplanes(expt))
end

function planelabel(expt::LPREExperiment, i)
    label = "τ = $(round(expt.recoverytimes[i]; sigdigits=3)) s"
    length(unique(filter(>(0), expt.concentrations))) > 1 || return label
    return label * ", c = $(expt.concentrations[i])"
end

function summarypanels(expt::LPREExperiment, peaks)
    return [ratepanel(expt, peaks, :Gamma2, "Γ₂"), ratepanel(expt, peaks, :Gamma1, "Γ₁"),
            intensitypanel(expt, peaks), transversepanel(expt, peaks),
            parampanel(peaks, :R1, "Diamagnetic R₁ / s⁻¹")]
end

# ---- the peak plot: height against recovery time -------------------------------------

crosssectionview(::LPREExperiment) = false
haslegend(::LPREExperiment) = true
legendposition(::LPREExperiment) = :rb
curvelabels(::LPREExperiment) = (xlabel="Recovery time / s", ylabel="Peak height")

function drawcurve!(ax, peak, expt::LPREExperiment)
    τ = expt.recoverytimes
    y = collect(peak.parameters[:amp].value[])
    e = collect(peak.parameters[:amp].uncertainty[])
    skip = skipset(expt)
    series = seriesplanes(expt)
    φp = fitvalues(peak, expt)
    A = peak.parameters[:amp0].value[][1]
    τs = range(0, 1.05 * maximum(τ), 100)
    for (j, (_, planes)) in enumerate(series)
        colour = seriescolour(j, length(series))
        used = filter(i -> i ∉ skip, planes)
        skipped = filter(i -> i ∈ skip, planes)
        errorbars!(ax, τ[used], y[used], e[used]; whiskerwidth=6, color=colour)
        scatter!(ax, τ[used], y[used]; color=colour)
        isempty(skipped) ||
            scatter!(ax, τ[skipped], y[skipped]; color=:transparent, strokecolor=:gray60,
                     strokewidth=1.5)
        peak.postfitted[] || continue
        c = expt.concentrations[first(planes)]
        lines!(ax, τs,
               [A * peakmodel(expt, (c=c, τ=t), φp[3], φp[4], φp[5:end])[3] for t in τs];
               color=colour)
    end
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::LPREExperiment)
    samples = expt.paired ? "diamagnetic, paramagnetic" :
              "concentrations $(join(expt.concentrations, ", "))"
    τ = join(round.(expt.recoverytimes; sigdigits=4), ", ")
    return "Analysis type: longitudinal PRE experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Samples: $samples\n" *
           "Recovery times / s: $τ\n" *
           "Experiment type: $(expt.expttype == :hsqc ? "HSQC" : "HMQC")\n" *
           "Transfer time Trelax / ms: $(round(1000 * expt.Trelax; digits=2))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end
