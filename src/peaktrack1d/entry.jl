# The entry point: resolve what was passed into spectra, then open the window.

"""
    peaktrack1d(expts; process=true, lb=1.0, zerofill=2, solvent=:auto, phase=:auto,
                labels=String[], reference="reference", refshift=0.0, coordinates=nothing,
                peaklist=nothing, radius=8, snr=5, gain=2, colormap=nothing,
                size=(1200, 850)) -> PeakTracker
    peaktrack1d(spectra; kwargs...)
    peaktrack1d(fids; lb=1.0, zerofill=2, solvent=:auto, phase=:auto, kwargs...)
    peaktrack1d(previous::PeakTracker; kwargs...)

Follow the positions of peaks across a series of 1D spectra, drawn as a stack. Each peak is
a multiplet of one or more lines: pick each line by painting along its trajectory, then
group the lines of a multiplet together. A multiplet's shift in each plane is the mean of
its lines' shifts. One multiplet, shown first, is the chemical shift reference: once it is
picked in every plane, the spectra can be aligned to it, and every shift is also reported
referenced to it.

The window blocks until it is closed, and the analysis is returned. Show it for a table of
multiplet centres, copy every line's shifts with `clipboard(p)`, get every pick with `planeresults(p)`, or
reopen the window with `peaktrack1d(p)`. **Save** writes the output folder named under the
plot, and **Load** restores the picks saved there.

# Arguments
- `expts`: an experiment folder or number, or a list or range of them (e.g. `11:42`),
  relative to the working directory. Each may be 1D or pseudo-2D.
- `spectra`: processed spectra, as one `NMRData` or a vector of them, sharing one chemical
  shift axis.
- `fids`: time-domain data from `loadfid`, as one `NMRData` or a vector of them.
- `process`: for `expts`, process the raw data (`true`, the default) or load the processed
  spectra in `pdata/1` (`false`).
- `lb`, `zerofill`, `solvent`: exponential line broadening (Hz), zero-filling factor, and
  solvent suppression as for NMRTools' `autoprocess`, when processing. The line broadening
  can be changed afterwards with the slider under the plot.
- `phase`: when processing, `:auto` to phase every experiment automatically, `:pdata` to
  take PHC0 and PHC1 from its processed data, a vector of `(p0, p1)` pairs in degrees, or
  the path of a `processing.csv` saved earlier. In a series of three or more, a first-order
  phase far from the series median is set to the median and the zero-order phase fitted
  again, and a zero-order phase far from the rest of the series takes its neighbours'.
- `labels`: names of multiplets to create before picking; more are created as you pick.
- `reference`: the name of the reference multiplet, picked in this analysis. Or, for ¹⁹F,
  ³¹P or another nucleus, a ¹H `PeakTracker` of the same samples with DSS picked: the
  shifts are then referenced indirectly to it on the IUPAC unified scale, using NMRTools'
  Ξ ratios, matching spectra one to one in order (see [`IndirectReference`](@ref)). The
  two must have the same number of spectra, recorded at the same field.
- `refshift`: the shift (ppm) at which a picked reference is placed when aligning and
  referencing.
- `coordinates`: what distinguishes the planes, one value per plane: a vector, a pair such
  as `"pH" => pH`, or a `NamedTuple` of vectors such as `(; pH, T)`. Written beside every
  shift in the output.
- `peaklist`: a `peaklist.csv` saved earlier, or an earlier `PeakTracker`, whose picks to
  start from. Each pick is placed exactly at its saved shift.
- `radius`: the search radius (pixels) for snapping and for the brush.
- `snr`: the smallest signal-to-noise ratio of a peak found by the brush. Picks placed by
  clicking are always accepted.
- `gain`: the initial vertical scale; `colormap` colours the traces by plane.

# Mouse
Only the current line responds; change it with the keys.
- **drag**: paint a trajectory for the current line. In every plane the brush crosses, the
  tallest peak whose maximum lies inside the shaded area, within that plane's strip of the
  stack, is picked. Narrow the brush with `[` and `]` to keep a strong neighbour out.
- **click**: pick the current line in the plane under the cursor.
- **drag a marker**: move it along its plane. **Right-click a marker**: delete it.
- **shift-drag**: erase the current line's picks under the stroke.
- **scroll**: zoom chemical shift. **Right-drag**: pan. **Ctrl-click**: reset the view.

Clicked and dragged picks snap to the tallest point within the search radius; `m` turns
snapping off (for overlapped peaks), and holding shift does the opposite of the setting.

# Keys
- `←` `→`: step through the lines, multiplet by multiplet, to the vacant slot at the end,
  where picking creates a new multiplet. `0`-`9`: select a multiplet; `0` is the reference,
  which can't be selected when it comes from a ¹H analysis.
- `r`: rename the current multiplet (type, then enter; escape cancels).
- `g`, then click a marker: join that line's multiplet to the current one. `u`: move the
  current line out into a multiplet of its own. `d`: delete the current line.
- `a`: align the spectra to the reference, once it is picked in every plane. This changes
  only the display; shifts are always stored as measured.
- `↑` `↓`: vertical scale. `[` `]`: search radius. `m`: snapping. `backspace`: undo.
- `s`: save the analysis to the output folder. `f`: save a figure of the region shown
  there, saving the analysis first if it hasn't been saved there from this window.
  `q`: quit.

# Examples
```julia
p = peaktrack1d(11:42)
p = peaktrack1d(11:42; process=false)               # the processed spectra instead
p = peaktrack1d(11:42; peaklist="out/peaklist.csv") # carry on from a saved analysis
p = peaktrack1d(p)                                  # reopen

h = peaktrack1d(11:42)                              # ¹H, with DSS picked as the reference
f = peaktrack1d(51:82; reference=h)                 # ¹⁹F of the same samples, on DSS
clipboard(p)
```

See also [`peaktrack2d`](@ref NMRAnalysis.GUI2D.peaktrack2d) for 2D spectra.
"""
function peaktrack1d end

function peaktrack1d(expts::AbstractVector{<:Union{Integer,AbstractString}};
                     process::Bool=true, zerofill=2, solvent=:auto, coordinates=nothing,
                     reference=DEFAULTREFERENCE, refshift=0.0, kwargs...)
    isempty(expts) && throw(ArgumentError("no experiments to track peaks in"))
    call = analysiscall("peaktrack1d", expts; process, zerofill, solvent,
                        coordinates=callcoordinates(coordinates), reference,
                        refshift)
    paths = string.(expts)
    data = process ? loadfid.(exptfolder.(paths)) : loadnmr.(paths)
    return peaktrack1d(data; zerofill, solvent, coordinates, reference, refshift, call,
                       kwargs...)
end

peaktrack1d(expt::Union{Integer,AbstractString}; kwargs...) = peaktrack1d([expt]; kwargs...)
peaktrack1d(spec::NMRData; kwargs...) = peaktrack1d([spec]; kwargs...)

function peaktrack1d(data::AbstractVector{<:NMRData}; lb=1.0, zerofill=2, solvent=:auto,
                     phase=:auto, call=analysiscall("peaktrack1d", data), kwargs...)
    if all(d -> dims(d, 1) isa TimeDimension, data)
        processing = Processing(data; lb, zerofill, solvent, phase)
        s = SpectrumSeries(processing)
    elseif all(d -> any(x -> x isa FrequencyDimension, dims(d)), data)
        processing = nothing
        s = SpectrumSeries(data)
    else
        throw(ArgumentError("expected processed spectra or time-domain data from loadfid, " *
                            "not a mixture"))
    end
    return openwindow(s, processing; call, kwargs...)
end

function peaktrack1d(p::PeakTracker; kwargs...)
    return openwindow(p.series, p.processing; coordinates=p.coordinates, peaklist=p,
                      reference=(multiplets(tracks(p))[REF].label, p.reference),
                      gain=p.view.gain[], call=p.call, outputdir=p.outputdir[], kwargs...)
end

"The experiment folder holding `path`, which may name its processed data instead."
function exptfolder(path::AbstractString)
    parts = splitpath(path)
    length(parts) ≥ 3 && lowercase(parts[end - 1]) == "pdata" &&
        return joinpath(parts[1:(end - 2)]...)
    return path
end

"Coordinates for the record of the call: as a `NamedTuple`, which has a literal form."
callcoordinates(::Nothing) = nothing
callcoordinates(c) = coordinatetable(c, nothing)
