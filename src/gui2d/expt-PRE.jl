"""
    pre2d(inputfilenames, paramagnetic_concs, expttype, Trelax)

!!! warning "Under development"
    PRE fitting is under active development. Results should be validated carefully before
    use in publication.

Start an interactive GUI for analysing paramagnetic relaxation enhancement (PRE)
experiments. The fitted PRE rate Γ (s⁻¹ per unit concentration) accounts for both
linebroadening and intensity attenuation across the series of spectra.

# Arguments
- `inputfilenames`: Vector of path strings to processed Bruker data directories, one per
  concentration point.
- `paramagnetic_concs`: Vector of paramagnetic agent concentrations, one per spectrum.
  Use `0` for the diamagnetic reference. For protein PREs (diamagnetic/paramagnetic pair),
  use `[0, 1]`; for solvent PRE titrations, provide the actual concentrations in any
  consistent unit (e.g. mM).
- `expttype`: Either `:hsqc` or `:hmqc`. In an HMQC experiment, PRE-induced broadening
  also affects the indirect-dimension linewidth (R2y); in an HSQC it affects only R2x.
- `Trelax`: Total magnetisation transfer time (seconds) during which PRE-induced
  relaxation accumulates. This is pulse-sequence-specific (typically twice the INEPT delay
  for HSQC, or longer for HMQC).

# Examples
```julia
# Protein PRE: diamagnetic and paramagnetic samples
pre2d(
    ["diamagnetic/pdata/1", "paramagnetic/pdata/1"],
    [0, 1],
    :hsqc,
    0.010
)

# Solvent PRE titration
pre2d(
    ["0mM/pdata/1", "1mM/pdata/1", "5mM/pdata/1", "10mM/pdata/1"],
    [0.0, 1.0, 5.0, 10.0],
    :hmqc,
    0.0089
)
```
"""
function pre2d(inputfilenames, paramagnetic_concs, expttype, Trelax)
    expt = PREExperiment(inputfilenames, paramagnetic_concs, expttype, Trelax)
    return gui!(expt)
end

"""
    PREExperiment

Paramagnetic relaxation enhancement experiment.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks  
- `paramagnetic_concs`: Vector of paramagnetic agent concentrations
- `expttype`: Either :hsqc or :hmqc
- `Trelax`: Relaxation time
"""
struct PREExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    paramagnetic_concs::Any
    expttype::Any
    Trelax::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function PREExperiment(specdata, peaks, paramagnetic_concs, expttype, Trelax)
        expt = new(specdata, peaks, paramagnetic_concs, expttype, Trelax,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}())
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

primaryparam(::PREExperiment) = :PRE

"""
    PREExperiment(inputfilenames, paramagnetic_concs, expttype, Trelax)

Create a PRE experiment from files and experimental parameters.
`expttype` should be `:hsqc` or `:hmqc`. `Trelax` is the timing during which relaxation
can occur during the sequence (magnetisation transfer delays etc) - this is specific to
the pulse sequence used.

Can be used to analyse solvent PREs, in which case concentrations should be specified -
or to analyse protein PREs, in which case concentrations should be set to 0 and 1 for
diamagnetic and paramagnetic states respectively.
"""
function PREExperiment(inputfilenames, paramagnetic_concs, exptexperimenttype, Trelax)
    exptexperimenttype in [:hsqc, :hmqc] ||
        throw(ArgumentError("Experiment type must be :hsqc or :hmqc"))

    specdata = preparespecdata(inputfilenames, paramagnetic_concs, PREExperiment)
    peaks = Observable(Vector{Peak}())

    return PREExperiment(specdata, peaks, paramagnetic_concs, exptexperimenttype, Trelax)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, paramagnetic_concs, ::Type{PREExperiment})
    @debug "Preparing spec data for relaxation experiment: $inputfilenames"

    spectra = map(loadnmr, inputfilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(paramagnetic_concs) do conc
        if conc == 0
            "Diamagnetic"
        else
            "Paramagnetic (conc = $conc)"
        end
    end

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

    # get initial values for amplitude
    x0, y0 = initialposition
    ix = findnearest(expt.specdata.x[1], x0)
    iy = findnearest(expt.specdata.y[1], y0)
    amp0 = expt.specdata.z[1][ix, iy]
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

"""Simulate single peak according to experiment type."""
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
        R2x = R2x0 + PRE * expt.paramagnetic_concs[i]
        if expt.expttype == :hmqc
            R2y = R2y0 + PRE * expt.paramagnetic_concs[i]
        else
            R2y = R2y0
        end
        amp = amp0 * exp(-PRE * expt.paramagnetic_concs[i] * expt.Trelax)

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
               "PRE: $(peak.parameters[:PRE].value[][1] ± peak.parameters[:PRE].uncertainty[][1]) s⁻¹ [conc⁻¹]\n" *
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
    return "Analysis type: PRE experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "PRE agent concs: $(join(expt.paramagnetic_concs, ", "))\n" *
           "Experiment type: $(expt.expttype==:hsqc ? "HSQC" : "HMQC")\n" *
           "Relaxation time: $(expt.Trelax)\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end
