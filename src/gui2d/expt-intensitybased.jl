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
    fit2d(inputfilenames)

Start an interactive GUI for peak analysis of a single 2D spectrum or a series of 2D
spectra. Each peak is fitted to a 2D Lorentzian lineshape; no physical model is applied
to the amplitudes across spectra.

Use this function to measure peak positions, linewidths, and amplitudes for downstream
analysis, or when none of the built-in physical models ([`relaxation2d`](@ref),
[`recovery2d`](@ref), [`modelfit2d`](@ref)) are appropriate.

# Arguments
- `inputfilenames`: A single path string or vector of path strings pointing to processed
  Bruker data directories (e.g. `"expno/pdata/1"`).

# Example
```julia
# Single spectrum
fit2d("109/pdata/1")

# Series of spectra (e.g. titration or temperature series)
fit2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"])
```
"""
function fit2d(inputfilenames)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    expt = IntensityExperiment(specdata,
                               peaks,
                               NoFitting())

    return gui!(expt)
end

"""
    relaxation2d(inputfilenames, relaxationtimes; skipplanes=nothing)

Start an interactive GUI for measuring R1 or R2 relaxation rates from a series of 2D
spectra. Peak amplitudes are fitted to a mono-exponential decay:

```math
I(\\tau) = A \\exp(-R\\tau)
```

where ``R`` is the relaxation rate (s⁻¹) and ``A`` is the peak amplitude. The software
does not distinguish R1 from R2 — the appropriate interpretation depends on the experiment.

# Arguments
- `inputfilenames`: Vector of path strings to processed Bruker data directories, one per
  relaxation delay.
- `relaxationtimes`: Vector of delay times in seconds, or a string giving a path to a
  text file containing the delays (one per line; lines beginning with `#` are ignored).

# Keyword Arguments
- `skipplanes`: Optional list of plane indices (1-based) to exclude from the exponential
  fit. All spectra are still loaded and displayed; skipped planes appear as open grey
  markers in the peak plot and are not used when fitting R or A. The full list of
  relaxation times must still be provided, including those for skipped planes.

# Example
```julia
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1"],
    [0.010, 0.030, 0.060, 0.100]
)

# Reading delays from a file
relaxation2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"], "vclist.txt")

# Omit the 3rd plane (e.g. corrupted or duplicate delay) from the fit
relaxation2d(
    ["11/pdata/1", "12/pdata/1", "13/pdata/1", "14/pdata/1"],
    [0.010, 0.030, 0.060, 0.100];
    skipplanes=[3]
)
```
"""
function relaxation2d(inputfilenames, relaxationtimes; skipplanes=nothing)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    tau = Float64[]
    if relaxationtimes isa String
        append!(tau, vec(readdlm(relaxationtimes; comments=true)))
    elseif relaxationtimes isa Vector
        for t in relaxationtimes
            if t isa String
                append!(tau, vec(readdlm(t; comments=true)))
            else
                append!(tau, t)
            end
        end
    end

    skip = isnothing(skipplanes) ? Int[] : collect(Int, skipplanes)
    if !isempty(skip)
        bad = filter(i -> i < 1 || i > length(tau), skip)
        isempty(bad) ||
            error("skipplanes indices out of range (got $bad for $(length(tau)) planes)")
    end

    expt = IntensityExperiment(specdata, peaks, ExponentialModel(), tau,
                               ModelFitVisualisation(); skipplanes=skip)

    return gui!(expt)
end

"""
    recovery2d(inputfilenames, relaxationtimes)

Start an interactive GUI for measuring longitudinal relaxation from an inversion recovery
or saturation recovery experiment. Peak amplitudes are fitted to a magnetisation recovery
model:

```math
I(\\tau) = A\\left(1 - C\\exp(-R\\tau)\\right)
```

where ``R`` is the recovery rate (s⁻¹), ``A`` is the equilibrium amplitude, and ``C``
is the recovery factor. For an ideal inversion recovery experiment ``C = 2``; for
saturation recovery ``C = 1``.

# Arguments
- `inputfilenames`: A single path string (pseudo-3D dataset) or vector of path strings
  (one file per delay) pointing to processed Bruker data directories.
- `relaxationtimes`: Vector of delay times in seconds, or a string giving a path to a
  text file containing the delays (one per line; lines beginning with `#` are ignored).

# Example
```julia
t = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
recovery2d("33/pdata/1", t)

# Reading delays from a file
recovery2d("33/pdata/1", "vdlist.txt")
```
"""
function recovery2d(inputfilenames, relaxationtimes)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    # First handle relaxation times
    tau = Float64[]
    if relaxationtimes isa String
        append!(tau, vec(readdlm(relaxationtimes; comments=true)))
    elseif relaxationtimes isa Vector
        for t in relaxationtimes
            if t isa String
                append!(tau, vec(readdlm(t; comments=true)))
            else
                append!(tau, t)
            end
        end
    end

    # specdata.zlabels .= map(t -> "τ = $t", tau)

    expt = IntensityExperiment(specdata,
                               peaks,
                               RecoveryModel(),
                               tau,
                               ModelFitVisualisation())

    return gui!(expt)
end

"""
    modelfit2d(inputfilenames, xvalues, equation, parameters)

Create an intensity analysis experiment with fitting to a custom equation.

# Arguments
- `inputfilenames`: String or vector of strings giving the input data files.
- `xvalues`: Vector of Float64 giving the x values for fitting, or string giving a filename
  from which to read the x values.
- `equation`: String giving the model equation to fit, e.g. `"A*sin(J*x)"`
- `parameters`: Vector of parameter name-value pairs giving initial parameter values,
  e.g. `["A"=>40., "J"=>0.5]`

# Example: J-modulation
```julia
modelfit2d(["112","113","114","115"],
    [0.1, 0.2, 0.3, 0.4],
    "A*sin(J*x)",
    ["A"=>40., "J"=>0.5])
```
"""
function modelfit2d(inputfilenames, xvalues, modelfunction::String,
                    parameters::Vector{Pair{String,Float64}}, xlabel="x")
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    # First handle relaxation times
    xval = Float64[]
    if xvalues isa String
        append!(xval, vec(readdlm(xvalues; comments=true)))
    elseif xvalues isa Vector
        for x in xvalues
            if x isa String
                append!(xval, vec(readdlm(x; comments=true)))
            else
                append!(xval, x)
            end
        end
    end

    # specdata.zlabels .= map(t -> "τ = $t", tau)
    model = CustomModel(modelfunction, parameters::Vector{Pair{String,Float64}}, xlabel)
    expt = IntensityExperiment(specdata,
                               peaks,
                               model,
                               xval,
                               ModelFitVisualisation())

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

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::IntensityExperiment, xbounds=nothing,
                   ybounds=nothing)
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]
    amp0 = peak.parameters[:amp].value[][1]

    for i in 1:nslices(expt)
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds[i]]
        y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds[i]]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]
        R2x = peak.parameters[:R2x].value[][i]
        R2y = peak.parameters[:R2y].value[][i]
        amp = peak.parameters[:amp].value[][i]

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

function postfit!(peak::Peak, expt::IntensityExperiment)
    return postfit!(peak, expt, expt.model)
end

function get_model_data(peak, expt::IntensityExperiment)
    return get_model_data(peak, expt, expt.model)
end

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
    skipped = idx in expt.skipplanes ? " [skipped]" : ""
    if length(expt.specdata.zlabels) == 1
        "Slice $idx of $(nslices(expt))$skipped"
    else
        "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))$skipped"
    end
end