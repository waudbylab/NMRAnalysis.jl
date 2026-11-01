"""
    pre2d(diamagnetic, paramagnetic; type, Trelax=:auto, peaklist=nothing) -> PREExperiment
    pre2d(inputfilenames; concentrations, type, Trelax=:auto, skipplanes=nothing,
          peaklist=nothing) -> PREExperiment

!!! warning "Under development"
    PRE fitting is under active development. Results should be validated carefully before
    use in publication.

Start an interactive GUI for analysing ¹H transverse paramagnetic relaxation enhancement
(PRE) experiments. The fitted PRE rate Γ₂ accounts for both line broadening and the
intensity lost during the coherence transfer delays.

Give a diamagnetic and a paramagnetic spectrum for a spin-labelled or metal-tagged protein,
and Γ₂ is in s⁻¹. Give a series of spectra with their `concentrations` for a solvent PRE
titration, and Γ₂ is in s⁻¹ per unit of concentration, in whatever unit you used.

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

# Solvent PRE titration, concentrations in mM
pre2d(["0mM/pdata/1", "1mM/pdata/1", "5mM/pdata/1", "10mM/pdata/1"];
      concentrations=[0, 1, 5, 10], type=:hmqc)
```

The window blocks until it is closed, and the analysis is returned (see [`results`](@ref)).
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
struct PREExperiment <: FixedPeakExperiment
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
        expttype in (:hsqc, :hmqc) ||
            throw(ArgumentError("type must be :hsqc or :hmqc, not $(repr(expttype))"))
        length(concentrations) == length(specdata.z) ||
            throw(ArgumentError("got $(length(concentrations)) concentrations for " *
                                "$(length(specdata.z)) spectra"))
        all(≥(0), concentrations) ||
            throw(ArgumentError("concentrations must not be negative"))
        used = [c for (i, c) in enumerate(concentrations) if i ∉ skipplanes]
        length(unique(used)) ≥ 2 ||
            throw(ArgumentError("PRE needs spectra at two or more concentrations"))
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

# One-bond couplings setting the nominal transfer delays: amide ¹H-¹⁵N and methyl ¹H-¹³C
const ONEBONDJ = Dict(N15 => 92.0, C13 => 125.0)

"""
    transfertime(Trelax, spec) -> Float64

`Trelax` in seconds, or for `:auto` the 1/¹J for the heteronucleus in F2 of `spec`: two
INEPTs of 1/(2J) in an HSQC, or the two 1/(2J) delays either side of t₁ in an HMQC.
"""
transfertime(Trelax::Real, spec) = Float64(Trelax)
function transfertime(Trelax::Symbol, spec)
    Trelax === :auto ||
        throw(ArgumentError("Trelax must be a time in seconds or :auto, " *
                            "not $(repr(Trelax))"))
    nuc = spec[F2Dim, :nucleus]
    haskey(ONEBONDJ, nuc) ||
        throw(ArgumentError("cannot set Trelax automatically for a $(label(spec, F2Dim)) " *
                            "indirect dimension: give it in seconds"))
    return 1 / ONEBONDJ[nuc]
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, zlabels, ::Type{PREExperiment})
    @debug "Preparing spec data for PRE experiment: $inputfilenames"

    spectra = map(loadnmr, inputfilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    return SpecData(spectra, x, y,
                    z ./ σ[1],
                    σ ./ σ[1],
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::PREExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)

    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(30.0)
    R2y0 = MaybeVector(15.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)

    # initial amplitude from the least paramagnetic spectrum
    i0 = argmin(expt.concentrations)
    x0, y0 = initialposition
    ix = findnearest(expt.specdata.x[i0], x0)
    iy = findnearest(expt.specdata.y[i0], y0)
    amp0 = expt.specdata.z[i0][ix, iy]
    amp0 = MaybeVector(amp0)
    amp = Parameter("Amplitude", amp0)

    # PRE
    Γ = Parameter("PRE", MaybeVector(10.0); minvalue=0.0, maxvalue=200.0)

    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp
    newpeak.parameters[:PRE] = Γ

    newpeak.postparameters[:PRE] = Parameter("PRE", 0.0)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""
Simulate single peak according to experiment type.

Γ₂ adds to the ¹H linewidth (and in an HMQC the multiple-quantum linewidth in F2) and
attenuates the amplitude by exp(-Γ₂ Trelax), which with the amplitude normalised to the
diamagnetic linewidths gives I_para/I_dia = R₂ exp(-Γ₂ t)/(R₂ + Γ₂) for peak heights
(Battiste & Wagner, Biochemistry 39, 5355 (2000); Iwahara, Tang & Clore, J. Magn. Reson.
184, 185 (2007)).
"""
function simulate!(z, peak::Peak, expt::PREExperiment, xbounds=nothing, ybounds=nothing)
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]
    amp0 = peak.parameters[:amp].value[][1]
    PRE = peak.parameters[:PRE].value[][1]

    for i in 1:nslices(expt)
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds[i]]
        y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds[i]]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]

        # apply PRE to linewidths and amplitude
        Γ = PRE * expt.concentrations[i]
        R2x = R2x0 + Γ
        R2y = expt.expttype == :hmqc ? R2y0 + Γ : R2y0
        amp = amp0 * exp(-Γ * expt.Trelax)

        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(2π * hz(x0, xaxis), R2x, 2π * hz(xs, xaxis),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x0 * R2y0) *
             NMRTools.NMRBase._lineshape(2π * hz(y0, yaxis), R2y, 2π * hz(ys, yaxis),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, ::PREExperiment)
    peak.postparameters[:PRE].uncertainty[] .= peak.parameters[:PRE].uncertainty[]
    peak.postparameters[:PRE].value[] .= peak.parameters[:PRE].value[]
    return peak.postfitted[] = true
end

# Γ₂ from a titration is per unit of a concentration whose unit was never given
paramunit(e::PREExperiment, name::Symbol) = name == :PRE && !e.paired ? "" : paramunit(name)

preunit(expt::PREExperiment) = expt.paired ? "s⁻¹" : "s⁻¹ per unit concentration"

"""Return descriptive text for slice idx."""
function slicelabel(expt::PREExperiment, idx)
    return "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::PREExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "PRE: $(peak.parameters[:PRE].value[][1] ± peak.parameters[:PRE].uncertainty[][1]) $(preunit(expt))\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "Amplitude: $(peak.parameters[:amp].value[][1] ± peak.parameters[:amp].uncertainty[][1])\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
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
