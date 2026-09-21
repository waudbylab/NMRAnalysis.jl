function componentnames(spec::NMRData)
    components = sample(spec, :sample, :components)
    isnothing(components) && return String[]
    return [component["name"] for component in components if haskey(component, "name")]
end
componentnames(expt::AbstractExperiment) = componentnames(expt.spec)

"""
    sampleconcentrations(spec::NMRData) -> Dict{String,Float64}

Extract a Dict mapping molecule names to concentrations from the NMR sample metadata.
"""
function sampleconcentrations(spec::NMRData)
    components = sample(spec, :sample, :components)
    isnothing(components) && return Dict{String,Float64}()
    result = Dict{String,Float64}()
    for component in components
        name = get(component, "name", nothing)
        conc = get(component, "concentration_or_amount", nothing)
        isnothing(name) && continue
        if isnothing(conc)
            @warn "Component \"$name\" has no concentration_or_amount defined — skipping"
            continue
        end
        result[name] = conc
    end
    return result
end
sampleconcentrations(expt::AbstractExperiment) = expt.sampleconcentrations

"""
    moleculeconcentration(model, expt::AbstractExperiment, role::Symbol) -> Float64

Look up the concentration of the molecule assigned to `role` (e.g. `:A`, `:X`)
for `expt`: first from the experiment's own sample metadata, falling back to
`model.concentrations` (populated interactively when metadata is missing).

Raises an informative `ArgumentError` naming the molecule, its role, and the
experiment when neither source has a value — e.g. because only some of the
loaded experiments were matched to a sample containing this molecule — rather
than letting a bare `KeyError` propagate from deep inside the fitting or
plotting code.
"""
function moleculeconcentration(model, expt::AbstractExperiment, role::Symbol)
    name = model.moleculemap[role]
    sc = sampleconcentrations(expt)
    haskey(sc, name) && return sc[name]
    haskey(model.concentrations, name) && return model.concentrations[name]
    return throw(ArgumentError("No concentration found for molecule \"$name\" (:$role) in " *
                               "experiment $(short_expt_path(expt)) — its sample metadata does " *
                               "not include \"$name\", and no fallback concentration was entered. " *
                               "Check that this experiment is matched to the correct sample."))
end

"""
    field_label(expt) -> Symbol

Convert a magnetic field strength in Tesla to a Symbol for use as a ComponentArray key.
E.g. `field_label(14.1)` → `:14p1T`.
"""
field_label(field_teslas::Float64) = Symbol(replace(string(field_teslas), "." => "p") * "T")
field_label(expt::AbstractExperiment) = field_label(expt.field_teslas)

"""
    resolvecalibration(given) -> (calibration, save)

Resolve the `calibration` argument an entry point was given, once per problem.

A ready-made [`B1Calibration`](@ref) is passed through and `nothing` stays `nothing`;
anything else is taken to be the nutation calibration experiment(s) to fit, which are
fitted here so that several experiments share one calibration rather than each refitting
it. The second return value writes that fit's own analysis into a folder when there is one
to write it into (see `calibrationanalysis`), and is `nothing` where nothing was fitted.
"""
resolvecalibration(::Nothing) = (nothing, nothing)
resolvecalibration(cal::B1Calibration) = (cal, nothing)
resolvecalibration(given) = calibrationanalysis(given)

"""
    b1calibration(given, spec, nuc) -> B1Calibration

The calibration one experiment should use: the one resolved for the problem, or the nominal
calibration read from that spectrum's own reference pulse, with an assumed B₁
inhomogeneity.

A calibration measured on another channel is a warning rather than an error, since
`exchange1d` may be given experiments on more than one nucleus.
"""
function b1calibration(given, spec, nuc)
    isnothing(given) && return B1Calibration(spec, nuc)
    calnuc = nucleus(given)
    isnothing(calnuc) || calnuc == nuc ||
        @warn "B₁ calibration was measured on $calnuc but is being applied to a " *
              "$nuc experiment"
    return given
end

include("expt-r1.jl")
include("expt-cest.jl")
include("expt-r1rho-onres.jl")
include("expt-r1rho-offres.jl")

"""
    load_experiment(filename) -> AbstractExperiment

Load an NMR experiment file and return the appropriate concrete experiment type
based on its annotations.

Dispatches on `annotations(spec, :experiment_type)` and `annotations(spec, :features)`:
- `"relaxation"` + `"R1"` → `R1Experiment`
- `"cest"` → `CESTExperiment`
- `"r1rho"` + `"on_resonance"`/`"off_resonance"` → `R1rhoOnResExperiment`/`R1rhoOffResExperiment`
"""
function load_experiment(filename; calibration=nothing)
    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("$filename has no annotations — cannot classify experiment"))

    types = annotations(spec, :experiment_type)
    features = annotations(spec, :features)

    if "relaxation" in types && "R1" in features
        return R1Experiment(filename)
    elseif "cest" in types
        return CESTExperiment(filename; calibration)
    elseif "r1rho" in types && "on_resonance" in features
        return R1rhoOnResExperiment(filename; calibration)
    elseif "r1rho" in types && "off_resonance" in features
        return R1rhoOffResExperiment(filename; calibration)
    else
        throw(ArgumentError("Cannot classify experiment $filename " *
                            "(types=$types, features=$features)"))
    end
end

"""
    ExchangeProblem(filenames::Vector{String}, model::AbstractModel; calibration=nothing)

Construct an ExchangeProblem by loading experiments from filenames.
Each file is classified and loaded via `load_experiment`.

`calibration` gives the B₁ field strengths and inhomogeneity every experiment is simulated
with: a [`B1Calibration`](@ref), or the nutation calibration experiment(s) to fit for one.
It is resolved here rather than per experiment, so calibration experiments are fitted once.
Without it each experiment falls back to its own reference pulse and an assumed 5%
inhomogeneity.
"""
function ExchangeProblem(filenames::Vector{String}, model::AbstractModel;
                         calibration=nothing)
    cal, savecalibration = resolvecalibration(calibration)
    experiments = AbstractExperiment[load_experiment(f; calibration=cal) for f in filenames]
    return ExchangeProblem(experiments, model, nothing, savecalibration)
end
