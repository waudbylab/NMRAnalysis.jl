"""
    PeakTrack1D

Interactive picking of peak positions across a series of 1D spectra, drawn as a stack. Each
entity is a multiplet made of one or more lines, and each line holds at most one pick per
plane. A reference multiplet, picked like any other, lets the spectra be aligned and the
shifts referenced.

The picking model (`tracks.jl`), the peak finding (`peaks.jl`) and the brush geometry
(`view.jl`) are plain functions of plain data, and are what the tests exercise. NMRData is
touched only when the spectra are loaded or processed (`series.jl`, `processing.jl`), and
Makie only in drawing and the window.
"""
module PeakTrack1D

using CairoMakie: CairoMakie
using Dates
using GLMakie
using InteractiveUtils
using NMRTools
using Statistics

import ..NMRAnalysis  # the module itself, for pkgversion(NMRAnalysis)
# shared output rules - see src/output.jl and docs/src/advanced/conventions.md
using ..NMRAnalysis: csvcolumn, csvvalue, safename, sanitizelabel, backupfolder, writetable,
                     shortpath
# the call that produced an analysis - see src/calls.jl
using ..NMRAnalysis: AnalysisCall, analysiscall, writereproduce, RECORD
using ..NMRAnalysis: windowname, outputfolder, capturewindow
# one table row per entity per plane, as for the 2D analyses
import ..GUI2D: planeresults

export peaktrack1d, PeakTracker, stackplot

include("series.jl")
include("tracks.jl")
include("referencing.jl")
include("peaks.jl")
include("processing.jl")
include("view.jl")
include("tracker.jl")
include("files.jl")
include("figures.jl")
include("gui.jl")
include("entry.jl")

end
