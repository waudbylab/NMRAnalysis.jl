"""
    IntensityExperiment <: FixedPeakExperiment

Generic experiment for measurement of intensity modulations across 2D spectra.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks  
"""
struct IntensityExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    x::Vector{Float64}
    model::FittingModel
    visualisation::VisualisationStrategy
    skipplanes::Vector{Int}

    function IntensityExperiment(specdata, peaks, model, xvalues=nothing,
                                 visualisation=CrossSectionVisualisation();
                                 skipplanes=Int[])
        if isnothing(xvalues)
            xvalues = 1.0 * collect(1:length(specdata.z))
        end
        expt = new(specdata, peaks,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}(),
                   xvalues, model, visualisation, collect(Int, skipplanes))
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        return expt
    end
end

visualisationtype(expt::IntensityExperiment) = expt.visualisation

# Primary derived parameter: amplitude when no model is fitted, the relaxation
# rate :R for exponential/recovery fits, otherwise the first model parameter.
function primaryparam(expt::IntensityExperiment)
    expt.model isa NoFitting && return :amp
    expt.model isa MethylCCRModel && return :S2tc  # derived order parameter × τc
    names = Symbol.(expt.model.param_names)
    return :R in names ? :R : first(names)
end

"""
    fit2d(inputfilenames; skipplanes=nothing) -> IntensityExperiment

Start an interactive GUI for peak analysis of a single 2D spectrum or a series of 2D
spectra. Each peak is fitted to a 2D Lorentzian lineshape; no physical model is applied
to the amplitudes across spectra.

Use this function to measure peak positions, linewidths, and amplitudes for downstream
analysis, or when none of the built-in physical models ([`relaxation2d`](@ref),
[`recovery2d`](@ref), [`modelfit2d`](@ref)) are appropriate. The window blocks until it
is closed, and the analysis is returned (see [`results`](@ref)).

# Arguments
- `inputfilenames`: A single path string or vector of path strings pointing to processed
  Bruker data directories (e.g. `"expno/pdata/1"`). Bruker experiment numbers work too.
- `skipplanes`: Planes (1-based) left out of the lineshape fit. They are still displayed,
  and their amplitudes still measured.

# Example
```julia
# Single spectrum
fit2d("109/pdata/1")

# Series of spectra (e.g. titration or temperature series)
fit2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"])
```
"""
function fit2d(inputfilenames; skipplanes=nothing)
    specdata = preparespecdata(asexptpath(inputfilenames), IntensityExperiment)
    skip = checkskipplanes(skipplanes, length(specdata.z))
    expt = IntensityExperiment(specdata, Observable(Vector{Peak}()), NoFitting();
                               skipplanes=skip)
    return gui!(expt)
end

"""
    relaxation2d(inputfilenames; relaxationtimes, ncyc, cycletime, skipplanes, prompt)
    relaxation2d(inputfilenames, relaxationtimes; kwargs...)

Start an interactive GUI for measuring R1 or R2 relaxation rates from a series of 2D
spectra. Peak amplitudes are fitted to a mono-exponential decay:

```math
I(\\tau) = A \\exp(-R\\tau)
```

where ``R`` is the relaxation rate (s⁻¹) and ``A`` is the peak amplitude. The software
does not distinguish R1 from R2 — the appropriate interpretation depends on the experiment.
The window blocks until it is closed, and the analysis is returned (see [`results`](@ref)).

The delays are taken from `relaxationtimes` if given; else from the loop counts `ncyc`
times the loop duration `cycletime`; else from the `relaxation.duration` annotation; else
from the `vdlist`; else from the `vclist` times `cycletime`; and if none of those is
available you are asked for them, as you are for a `cycletime` that is needed and not given.

# Arguments
- `inputfilenames`: A single path string (pseudo-3D dataset) or vector of path strings
  (one file per delay) pointing to processed Bruker data directories. Bruker experiment
  numbers work too.

# Keyword Arguments
- `relaxationtimes`: Delays in seconds, one per plane, or the path of a text file holding
  one per line (lines beginning with `#` are ignored).
- `ncyc`: Loop counts, one per plane, or the path of a file holding them, for an experiment
  whose delay is a number of loops of fixed duration.
- `cycletime`: The duration in seconds of one loop, which multiplies `ncyc` or the `vclist`.
- `skipplanes`: Planes (1-based) to exclude from the fitting. All spectra are still loaded
  and displayed; skipped planes appear as open grey markers in the peak plot and are not
  used when fitting R or A. The full list of relaxation times must still be provided,
  including those for skipped planes.
- `prompt`: Whether to ask for anything that cannot be found (default: when interactive).

# Example
```julia
relaxation2d(["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1"];
             relaxationtimes=[0.010, 0.030, 0.060, 0.100])

# A pseudo-3D dataset whose delays are in its vdlist or annotations
relaxation2d("11/pdata/1")

# A CPMG-type R2 experiment counting loops in a vclist, each 16 ms long
relaxation2d("12/pdata/1"; cycletime=0.016)

# Omit the 3rd plane (e.g. corrupted or duplicate delay) from the fit
relaxation2d("11/pdata/1"; skipplanes=[3])
```
"""
function relaxation2d(inputfilenames; relaxationtimes=nothing, ncyc=nothing,
                      cycletime=nothing, skipplanes=nothing, prompt::Bool=isinteractive())
    specdata = preparespecdata(asexptpath(inputfilenames), IntensityExperiment)
    tau = relaxationdelays(specdata; relaxationtimes, ncyc, cycletime, prompt)
    skip = checkskipplanes(skipplanes, length(tau))
    model = ExponentialModel()
    expt = IntensityExperiment(specdata, Observable(Vector{Peak}()), model, tau,
                               ModelFitVisualisation(); skipplanes=skip)
    return gui!(expt)
end

function relaxation2d(inputfilenames, relaxationtimes; kwargs...)
    return relaxation2d(inputfilenames; relaxationtimes, kwargs...)
end

"""
    recovery2d(inputfilenames; relaxationtimes, ncyc, cycletime, skipplanes, prompt)
    recovery2d(inputfilenames, relaxationtimes; kwargs...)

Start an interactive GUI for measuring longitudinal relaxation from an inversion recovery
or saturation recovery experiment. Peak amplitudes are fitted to a magnetisation recovery
model:

```math
I(\\tau) = A\\left(1 - C\\exp(-R\\tau)\\right)
```

where ``R`` is the recovery rate (s⁻¹), ``A`` is the equilibrium amplitude, and ``C``
is the recovery factor. For an ideal inversion recovery experiment ``C = 2``; for
saturation recovery ``C = 1``. The window blocks until it is closed, and the analysis is
returned (see [`results`](@ref)).

The delays are resolved as for [`relaxation2d`](@ref), whose keyword arguments this shares.

# Example
```julia
t = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
recovery2d("33/pdata/1"; relaxationtimes=t)

# Reading delays from a file
recovery2d("33/pdata/1"; relaxationtimes="vdlist.txt")
```
"""
function recovery2d(inputfilenames; relaxationtimes=nothing, ncyc=nothing,
                    cycletime=nothing, skipplanes=nothing, prompt::Bool=isinteractive())
    specdata = preparespecdata(asexptpath(inputfilenames), IntensityExperiment)
    tau = relaxationdelays(specdata; relaxationtimes, ncyc, cycletime, prompt)
    skip = checkskipplanes(skipplanes, length(tau))
    expt = IntensityExperiment(specdata, Observable(Vector{Peak}()), RecoveryModel(), tau,
                               ModelFitVisualisation(); skipplanes=skip)
    return gui!(expt)
end

function recovery2d(inputfilenames, relaxationtimes; kwargs...)
    return recovery2d(inputfilenames; relaxationtimes, kwargs...)
end

"""
    modelfit2d(inputfilenames, xvalues, equation, parameters, xlabel="x";
               skipplanes=nothing, prompt=isinteractive())

Create an intensity analysis experiment with fitting to a custom equation. The window
blocks until it is closed, and the analysis is returned (see [`results`](@ref)).

# Arguments
- `inputfilenames`: String or vector of strings giving the input data files.
- `xvalues`: Vector of Float64 giving the x values for fitting, or string giving a filename
  from which to read the x values. `nothing` asks for them.
- `equation`: String giving the model equation to fit, e.g. `"A*sin(J*x)"`
- `parameters`: Vector of parameter name-value pairs giving initial parameter values,
  e.g. `["A"=>40., "J"=>0.5]`
- `skipplanes`: Planes (1-based) to exclude from the fitting.

# Example: J-modulation
```julia
modelfit2d(["112","113","114","115"],
    [0.1, 0.2, 0.3, 0.4],
    "A*sin(J*x)",
    ["A"=>40., "J"=>0.5])
```
"""
function modelfit2d(inputfilenames, xvalues, modelfunction::String,
                    parameters::Vector{Pair{String,Float64}}, xlabel="x";
                    skipplanes=nothing, prompt::Bool=isinteractive())
    specdata = preparespecdata(asexptpath(inputfilenames), IntensityExperiment)
    n = length(specdata.z)
    xval = isnothing(xvalues) ? askvector("x values", n; prompt) : readvalues(xvalues)
    checklength(xval, n, "x values")
    skip = checkskipplanes(skipplanes, n)
    model = CustomModel(modelfunction, parameters, xlabel)
    expt = IntensityExperiment(specdata, Observable(Vector{Peak}()), model, xval,
                               ModelFitVisualisation(); skipplanes=skip)
    return gui!(expt)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, ::Type{IntensityExperiment})
    @debug "Preparing spec data for intensity experiment: $inputfilenames"

    spec, x, y, z, σ, zlabels = if inputfilenames isa String
        # load a single file
        spec, x, y, z, σ = loadspecdata(inputfilenames, IntensityExperiment)
        (SingleElementVector(spec),
         SingleElementVector(x),
         SingleElementVector(y),
         z ./ σ,
         SingleElementVector(1),
         SingleElementVector(choptitle(label(spec))))
    elseif inputfilenames isa Vector{String}
        # load multiple files
        tmp = loadspecdata.(inputfilenames, IntensityExperiment)
        spec = []
        x = []
        y = []
        z = []
        σ = []
        zlabels = []
        for t in tmp
            n = length(t[4])
            if n == 1 # z is a single slice
                push!(spec, t[1])
                push!(x, t[2])
                push!(y, t[3])
                push!(z, t[4][1])
                push!(σ, t[5])
                push!(zlabels, choptitle(label(t[1])))
            else
                append!(spec, fill(t[1], n))
                append!(x, fill(t[2], n))
                append!(y, fill(t[3], n))
                append!(z, t[4])
                append!(σ, fill(t[5], n))
                append!(zlabels, fill(choptitle(label(t[1])), n))
            end
        end
        map(MaybeVector, (spec, x, y, z ./ σ[1], σ, zlabels))
    end

    return SpecData(spec, x, y, z, σ, zlabels)
end

# load the NMR data and prepare the SpecData object
function loadspecdata(inputfilename, ::Type{IntensityExperiment})
    @debug "Loading spec data for intensity experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)

    dat = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    z = if ndims(spec) == 3
        eachslice(dat; dims=3)
    else
        [dat]
    end

    return spec, x, y, z, σ
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::IntensityExperiment, initialposition::Point2f, label="",
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
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)

    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    # Add post-parameters based on model type
    setup_post_parameters!(newpeak, expt.model)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

# No post-parameters needed for NoFitting
function setup_post_parameters!(::Peak, ::NoFitting) end

# Add post-parameters for parametric models
function setup_post_parameters!(peak::Peak, model::ParametricModel)
    for name in model.param_names
        peak.postparameters[Symbol(name)] = Parameter(name, 0.0)
    end
end

function postfit!(peak::Peak, expt::IntensityExperiment)
    return postfit!(peak, expt, expt.model)
end

function get_model_data(peak, expt::IntensityExperiment)
    return get_model_data(peak, expt, expt.model)
end

fittedamplitudes(peak, expt::IntensityExperiment) = fittedamplitudes(peak, expt, expt.model)

function peakinfotext(expt::IntensityExperiment, idx)
    if idx == 0
        return "No peak selected"
    end

    peak = expt.peaks[][idx]
    if !peak.postfitted[]
        return "Peak: $(peak.label[])\nNot fitted"
    end

    # Common peak information
    info = ["Peak: $(peak.label[])",
            "",
            "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm",
            "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm",
            "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹",
            "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"]

    # Add model-specific parameters
    append!(info, model_parameter_text(peak, expt.model))

    return join(info, "\n")
end

function experimentinfo(expt::IntensityExperiment)
    info = ["Analysis type: Intensity",
            "Model: $(typeof(expt.model))",
            "Filename: $(expt.specdata.nmrdata[1][:filename])",
            "Number of peaks: $(length(expt.peaks[]))",
            "Experiment title: $(expt.specdata.nmrdata[1][:title])"]

    append!(info, model_info_text(expt.model, expt.x))

    isempty(expt.skipplanes) ||
        push!(info, "Skipped planes: $(join(expt.skipplanes, ", "))")

    return join(info, "\n")
end

get_model_xlabel(expt::IntensityExperiment) = expt.model.xlabel
function get_model_ylabel(expt::IntensityExperiment)
    return expt.model isa MethylCCRModel ? "|Iₐ / I_b|" : "Peak amplitude"
end

function slicelabel(expt::IntensityExperiment, idx)
    if length(expt.specdata.zlabels) == 1
        "Slice $idx of $(nslices(expt))"
    else
        "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
    end
end
