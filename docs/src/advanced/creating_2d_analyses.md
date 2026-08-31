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

## Required Implementation

Concrete subtypes must implement both the core analysis functions below and the visualization functions documented separately. Note that many other functions have default implementations that can be used unless special behaviour is needed.

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

### `simulate!(z, peak, expt, [xbounds], [ybounds])`
Simulate a single peak in the experiment.

**Arguments:**
- `z`: Array to store simulation results
- `peak`: The peak to simulate
- `expt`: The experiment
- `xbounds`: Optional bounds for x dimension simulation
- `ybounds`: Optional bounds for y dimension simulation

## Default Implementations

The abstract type provides several default implementations that can be used as-is or overridden when needed:

### `mask!(z, peak, expt)`
Default implementation generates an elliptical mask for each peak using `maskellipse!`. Only needs to be overridden if a different peak shape is required.

**Arguments:**
- `z`: Array to store mask results
- `peak`: The peak to mask
- `expt`: The experiment

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
- `mask!(z, peaks, expt)`: Calculate peak masks and update internal specdata
- `simulate!(z, peaks, expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the peaks in the experiment
- `fit!(cluster, expt)`: Fit a specific cluster of peaks
- `checktouched!(expt)`: Check which clusters have been modified

### Observable Setup
- `setupexptobservables!(expt)`: Set up reactive behaviours for experiment observables

### Required Visualization Functions

Concrete subtypes must implement the following visualization functions:

### `makepeakplot!(gui, state, expt)`
Create the interactive peak plot in the GUI context. This function is crucial for real-time visualization and interaction.

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
- `bounds(mask)`: Calculate bounds from a mask
- `peak_plot_data(peak, expt)`: Extract plotting data for a single peak
- `plot_peak!(ax, peak, expt)`: Plot a single peak's data

## Implemented Experiment Types

The codebase includes implementations for several specific experiment types:
- `RelaxationExperiment`: For relaxation measurements
- `HetNOEExperiment`: For heteronuclear NOE measurements
- `PREExperiment`: For paramagnetic relaxation enhancement measurements

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
       
       # Experiment-specific fields
       my_special_parameter
   end
   ```

3. **Constructor**
   - Create a constructor that initializes all required fields
   - Set up observables using `setupexptobservables!`
   - Initialize experiment-specific parameters

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
   - Implement `simulate!` for your specific peak shapes/behavior
   - Implement `postfit!` to calculate final parameters from fit results

5. **Information Functions**
   - Implement `slicelabel` for spectrum navigation
   - Implement `peakinfotext` to show fit results
   - Implement `experimentinfo` to show experiment details

6. **Visualization Functions**
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
   - Override `mask!` only if you need non-elliptical peak shapes
   - Override clustering functions only if you need special peak grouping
   - Override other default implementations only if needed

Remember:
- Post-fit parameters (`postparameters`) are what get saved in results files
- Peak parameters (`parameters`) are used during fitting
- Use the existing implementations (RelaxationExperiment, HetNOEExperiment, PREExperiment) as templates
- Most functionality can be inherited from the abstract type