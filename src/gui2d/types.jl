"""
Abstract type representing an NMR experiment.

Subtypes must inherit from either FixedPeakExperiment or MovingPeakExperiment.
See guide for implementing new experiment types.
"""
abstract type Experiment end

abstract type FixedPeakExperiment <: Experiment end
abstract type MovingPeakExperiment <: Experiment end

abstract type FittingModel end
struct NoFitting <: FittingModel end
abstract type ParametricModel <: FittingModel end

struct SpecData
    nmrdata::Any
    x::Any
    y::Any
    z::Any
    σ::Any
    zlabels::Any
    zfit::Any
    mask::Any
end

struct Parameter
    label::Any
    value::Any
    uncertainty::Any
    initialvalue::Any
    minvalue::Any
    maxvalue::Any
end

struct Peak
    parameters::OrderedDict{Symbol,Parameter}
    label::Observable{String}
    touched::Observable{Bool}
    xradius::Observable{Float64}
    yradius::Observable{Float64}
    postparameters::OrderedDict{Symbol,Parameter}
    postfitted::Observable{Bool}
end

abstract type VisualisationStrategy end

struct CrossSectionVisualisation <: VisualisationStrategy end
struct ModelFitVisualisation <: VisualisationStrategy end
