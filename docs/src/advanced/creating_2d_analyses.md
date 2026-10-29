# 2D Experiments

An abstract type representing a type of experimental NMR analysis. The type system is split between fixed-peak experiments (where peak positions are constant between spectra) and peak-tracking experiments (where peak positions can vary).

## Required Fields

All concrete subtypes must include the following fields:
- `specdata`: A `SpecData` object containing the observed and simulated data plus mask
- `peaks`: An `Observable` list of peaks in the experiment
- `clusters`: An `Observable` list of clusters of peaks
- `touched`: An `Observable` list of touched clusters
- `isfitting`: An `Observable` boolean indicating if real-time fitting is active
- `xradius`: An `Observable` number defining peak detection radius in x dimension
- `yradius`: An `Observable` number defining peak detection radius in y dimension
- `state`: An `Observable` dictionary of GUI state variables
- `skipplanes`: A `Vector{Int}` of planes left out of the fitting (see `skipset`)

## Required Implementation

Concrete subtypes must implement both the core analysis functions below and the visualisation functions documented separately. Note that many other functions have default implementations that can be used unless special behaviour is needed.

### Type Hierarchy

The `Experiment` type has two immediate subtypes that handle position behaviour:
- `FixedPeakExperiment`: Implementation where `hasfixedpositions(expt) = true`
- `MovingPeakExperiment`: Implementation where `hasfixedpositions(expt) = false`

Concrete experiment types should inherit from one of these intermediate types rather than directly from `Experiment`.

### `addpeak!(expt, position, [label], [xradius], [yradius])`
Add a new peak to the experiment at the specified position.

**Arguments:**
- `expt`: The experiment
- `position`: A `Point2f` specifying the (x,y) coordinates of the peak
- `label`: Optional string label for the peak (defaults to auto-generated)
- `xradius`: Optional peak radius in x dimension (defaults to experiment default)
- `yradius`: Optional peak radius in y dimension (defaults to experiment default)

Every peak needs the parameters `:x`, `:y`, `:R2x`, `:R2y` and `:amp`, which the shared
lineshape fit in `fitting.jl` works on. A fixed-peak experiment holds one position and pair
of linewidths for the whole series and one amplitude per plane; a moving-peak experiment
holds all five per plane.

## Lineshape fitting

The fitting is shared by every experiment and needs nothing from a new one. Each peak is
the product of two lineshapes, scaled so that its amplitude is the peak height, and the
model is linear in the amplitudes. They are solved for exactly by least squares at every
step of the nonlinear fit (variable projection), so for a fixed-peak experiment only four
parameters per peak are optimised, however many planes there are; a moving-peak experiment
is fitted plane by plane. Planes in `skipplanes` do not constrain the shapes of a
fixed-peak experiment, but their amplitudes are still measured.

`fit!(expt)` fits the touched clusters in a background task, in parallel when Julia has
threads. Nothing is written into a peak until its cluster's fit is complete, so a fit that
is superseded or cancelled changes nothing. Each peak records how its last fit ended, in
`fitstatus(peak)`: `:converged`, `:maxiter`, `:bound` or `:timeout`.

## Default Implementations

The abstract type provides several default implementations that can be used as-is or overridden when needed:

### `maskplane!(m, peak, expt, i)`
Default implementation marks an elliptical fitting region for the peak in plane `i` using `maskellipse!`. Only needs to be overridden if a different region shape is required.

**Arguments:**
- `m`: The plane's mask, a `BitMatrix`
- `peak`: The peak to mask
- `expt`: The experiment
- `i`: The plane

### `postfit!(peak, expt)`
Perform additional fitting operations after spectrum fitting.

**Arguments:**
- `peak`: The peak to post-fit
- `expt`: The experiment

### `slicelabel(expt, idx)`
Generate a label for a specific slice/plane of the experiment.

**Arguments:**
- `expt`: The experiment
- `idx`: The slice index

**Returns:**
- A string label for the slice

### `peakinfotext(expt, idx)`
Generate information text about a specific peak.

**Arguments:**
- `expt`: The experiment
- `idx`: The peak index

**Returns:**
- A string containing formatted peak information

### `experimentinfo(expt)`
Generate information text about the experiment.

**Returns:**
- A string containing formatted experiment information

### `completestate!(state, expt)`
Set up observables for the GUI state.

**Arguments:**
- `state`: The GUI state dictionary
- `expt`: The experiment

## Functions Handled by Abstract Type

The following functions are implemented generically and do not need to be reimplemented by concrete subtypes:

### Peak Management
- `nslices(expt)`: Get the number of slices in the experiment
- `npeaks(expt)`: Get the number of peaks in the experiment
- `movepeak!(expt, idx, newpos)`: Move a peak to a new position
- `deletepeak!(expt, idx)`: Delete a specific peak
- `deleteallpeaks!(expt)`: Delete all peaks

### Data Processing
- `mask!(expt)`: Calculate peak masks and update internal specdata
- `simulate!(z, peaks, expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the touched clusters of peaks in the experiment
- `fitcluster!(peaks, expt, check)`: Fit one cluster of overlapping peaks
- `continuefit!(expt)`: Continue the fits stopped at the time or iteration limit
- `simulate!(z, peak, expt)`: Add a peak's fitted lineshape to every plane of `z`
- `batchupdate(f, expt)`: Make several changes to the peaks, then mask, cluster and fit once
- `checktouched!(expt)`: Check which clusters have been modified

### Observable Setup
- `setupexptobservables!(expt)`: Set up reactive behaviours for experiment observables

### Required Visualisation Functions

Concrete subtypes must implement the following visualisation functions:

### `makepeakplot!(gui, state, expt)`
Create the interactive peak plot in the GUI context. This function is crucial for real-time visualisation and interaction.

**Arguments:**
- `gui`: The GUI context containing plot panels
- `state`: The GUI state dictionary
- `expt`: The experiment

### `save_peak_plots!(expt, folder)`
Save publication-quality plots for all peaks to a specified folder.

**Arguments:**
- `expt`: The experiment
- `folder`: String path to the output folder

### Utility Functions
- `get_model_data(peak, expt)`: The observed and fitted series for a single peak
- `plot_peak!(panel, peak, expt)`: Plot a single peak's data

## Implemented Experiment Types

The codebase includes implementations for several specific experiment types:
- `RelaxationExperiment`: For relaxation measurements
- `HetNOEExperiment`: For heteronuclear NOE measurements

Each implementation specialises the simulation and fitting behaviour for its specific experiment type while inheriting the common functionality from the abstract type.

## Guide: Creating a New Experiment Type

Here's a step-by-step guide to implementing a new type of NMR experiment:

1. **Choose Base Type**
   - Inherit from `FixedPeakExperiment` if peak positions are constant between spectra
   - Inherit from `MovingPeakExperiment` if peak positions can vary

2. **Define Structure**
   ```julia
   struct MyNewExperiment <: FixedPeakExperiment
       # Required fields
       specdata
       peaks
       clusters
       touched
       isfitting
       xradius
       yradius
       state
       skipplanes::Vector{Int}

       # Experiment-specific fields
       my_special_parameter
   end
   ```

3. **Constructor**
   - Create a constructor that initialises all required fields
   - Set up observables using `setupexptobservables!`
   - Initialise experiment-specific parameters

4. **Core Analysis Functions**
   - Implement `addpeak!` to set up experiment-specific peak parameters
     ```julia
     function addpeak!(expt::MyNewExperiment, initialposition::Point2f, label="")
         # Create basic peak
         newpeak = Peak(initialposition, label)
         
         # Add experiment-specific parameters
         newpeak.parameters[:my_param] = Parameter("My Parameter", initial_value)
         
         # Add post-fit parameters that will be saved in results
         newpeak.postparameters[:final_result] = Parameter("Final Result", 0.0)
         
         push!(expt.peaks[], newpeak)
         notify(expt.peaks)
     end
     ```
   - Implement `postfit!` to calculate final parameters from fit results, leaving out
     the planes in `skipset(expt)`

5. **Information Functions**
   - Implement `slicelabel` for spectrum navigation
   - Implement `peakinfotext` to show fit results
   - Implement `experimentinfo` to show experiment details

6. **Visualisation Functions**
   - Implement `makepeakplot!` for the interactive GUI
     ```julia
     function makepeakplot!(gui, state, expt::MyNewExperiment)
         # Create appropriate plot(s) for your data
         gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1,1],
                                    xlabel="My X Label",
                                    ylabel="My Y Label")
         # Add plot elements
         plot!(ax, ...)
     end
     ```
   - Implement `save_peak_plots!` for publication figures
     ```julia
     function save_peak_plots!(expt::MyNewExperiment, folder::AbstractString)
         CairoMakie.activate!()
         for peak in expt.peaks[]
             fig = Figure()
             # Create publication-quality figure
             save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
         end
         GLMakie.activate!()
     end
     ```

7. **Optional Overrides**
   - Override `maskplane!` only if you need non-elliptical fitting regions
   - Override clustering functions only if you need special peak grouping
   - Override other default implementations only if needed

Remember:
- Post-fit parameters (`postparameters`) are what get saved in results files
- Peak parameters (`parameters`) are used during fitting
- Use the existing implementations (IntensityExperiment, HetNOEExperiment) as templates
- Resolve each parameter the entry point needs as the 1D analyses do: an explicit keyword,
  then an annotation, then an acquisition parameter, then `ask` (see `src/prompts.jl`)
- The entry point returns `gui!(expt)`, which blocks until the window closes and returns
  the experiment
- Most functionality can be inherited from the abstract type