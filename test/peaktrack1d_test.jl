# PeakTrack1D: the picking model, peak finding, brush geometry and output tables, with no
# window and no data files. Spectra are built in memory as a `SpectrumSeries`, the seam
# below which nothing touches NMRData. Processing and phasing need real FIDs and are left to
# tests alongside example data.

using NMRAnalysis
using NMRAnalysis.PeakTrack1D: SpectrumSeries, Pick, Multiplet, Tracks, REF, StackView,
                               PeakTracker
using NMRAnalysis.PeakTrack1D: climb, vertex, snap, exact, snr, noiselevel, displayscale,
                               sameaxis, planesfrom
using NMRAnalysis.PeakTrack1D: cleanlabel, multipletof, multipletlabels, linesof, orderedlines, navorder,
                               centre, linecentre, referenceshifts, vacantname, addmultiplet!,
                               setpick!, erase!, group!, ungroup!, deleteline!, relabel!,
                               haspicks, picksof, npicks, pickat, tracks
using NMRAnalysis.PeakTrack1D: densify, painted, track, overlay, hit, nearestrow, location
using NMRAnalysis.PeakTrack1D: coordinatetable, seriestable, peaklisttable, clipboardtext,
                               isreferenced, widetable,
                               linetable, referencedshift,
                               restore, readtable, saveanalysis, exptfolder, callcoordinates,
                               p1outliers, P1TOLERANCE
using NMRAnalysis: writetable
using CairoMakie: Figure, Axis, Observable, Point2f, Vec2f
using Test

# ---- synthetic spectra --------------------------------------------------------

const δaxis = collect(range(10.0, 0.0; length=2001))   # 0.005 ppm per point, high to low

lorentz(a, δ0; hw=0.01) = @. a * hw^2 / (hw^2 + (δaxis - δ0)^2)

# a reference that drifts, a singlet that moves downfield, and a doublet (lines 0.04 ppm
# apart) that moves upfield
dss(i) = 0.5 + 0.002i
singlet(i) = 3.0 + 0.05(i - 1)
doublet(i) = 7.0 - 0.03(i - 1)

function makeseries(n=4)
    Y = reduce(hcat,
               [lorentz(50, dss(i)) .+ lorentz(100, singlet(i)) .+
                lorentz(80, doublet(i) + 0.02) .+ lorentz(80, doublet(i) - 0.02) for i in 1:n])
    return SpectrumSeries(δaxis, Y, ones(n), fill("synthetic", n), 1:n; nucleus="1H",
                          title="test")
end

"Tracks with the reference, the singlet and both doublet lines snapped in every plane."
function picked(s; group=true)
    t = Tracks(["DSS"])
    for i in eachindex(s)
        setpick!(t, 1, snap(s, i, dss(i), 0.01))
    end
    a = addmultiplet!(t, "MES")
    b = addmultiplet!(t, "MOPS")
    c = addmultiplet!(t, "peak3")
    for i in eachindex(s)
        setpick!(t, a, snap(s, i, singlet(i), 0.01))
        setpick!(t, b, snap(s, i, doublet(i) + 0.02, 0.01))
        setpick!(t, c, snap(s, i, doublet(i) - 0.02, 0.01))
    end
    group && group!(t, b, c)
    return t
end

function tracker(s, t; coordinates="pH" => [6.0, 6.5, 7.0, 7.5][1:length(s)], refshift=0.0)
    fig = Figure()
    ax = Axis(fig[1, 1])
    return PeakTracker(s, coordinatetable(coordinates, length(s)), Observable(t), nothing,
                       refshift, Observable(false), StackView(s), fig, ax, nothing,
                       Observable("out"))
end

# ---- the data -----------------------------------------------------------------

@testset "series" begin
    s = makeseries()
    @test length(s) == 4
    @test collect(eachindex(s)) == 1:4
    @test planesfrom(s, 3) == [3]
    @test noiselevel(s, 2) == 1.0
    @test noiselevel([1.0, -1.0, 1.0, -1.0]) ≈ 1.4826
    @test sameaxis(δaxis, copy(δaxis))
    @test !sameaxis(δaxis, δaxis .+ 0.001)
    @test !sameaxis(δaxis, δaxis[1:(end - 1)])
    # the 90th percentile of the peak heights: between the doublet lines and the singlet,
    # and nowhere near a single dominant line
    @test 50 < displayscale(s) ≤ 101
    @test_throws DimensionMismatch SpectrumSeries(δaxis, zeros(10, 2), ones(2), ["a", "b"],
                                                  [1, 1])
end

@testset "peak finding" begin
    @test climb([0, 1, 3, 2, 0], 1) == 3
    @test climb([0, 1, 3, 2, 0], 5) == 3
    x = [1.0, 2.0, 3.0]
    y = @. 5 - (x - 2.3)^2
    @test all(vertex(x, y, 2) .≈ (2.3, 5.0))
    @test vertex(x, y, 1) == (1.0, y[1])   # no neighbour to interpolate with at an end

    s = makeseries()
    q = snap(s, 2, singlet(2) + 0.008, 0.01)
    @test q.plane == 2
    @test q.shift ≈ singlet(2) atol = 1e-3
    @test q.height ≈ 100 atol = 1
    @test isnothing(snap(s, 2, 20.0, 0.01))   # nothing within the window
    q = exact(s, 1, 3.0123)
    @test q.shift == 3.0123
    @test q.height == s.intensities[argmin(abs.(δaxis .- 3.0123)), 1]
    @test snr(s, Pick(1, 3.0, 42.0)) == 42.0
end

# ---- the picking model ----------------------------------------------------------

@testset "tracks" begin
    t = Tracks(["DSS", "MES"])
    @test multipletlabels(t) == ["DSS", "MES"]
    @test linesof(t, 2) == [2]
    @test t.nextline == 3
    @test vacantname(multipletlabels(t)) == "peak1"
    @test vacantname(["peak1", "peak3"]) == "peak2"
    @test cleanlabel(" a_b,c ") == "a-bc"
    l = addmultiplet!(t, "peak1")
    @test l == 3 && multipletof(t, 3) == 3

    setpick!(t, 2, Pick(1, 3.0, 10.0))
    setpick!(t, 2, Pick(1, 3.1, 11.0))   # replaces the earlier pick in that plane
    setpick!(t, 2, nothing)              # a failed snap places nothing
    @test npicks(t) == 1
    @test pickat(t, 2, 1).shift == 3.1
    @test haspicks(t, 2) && !haspicks(t, 3)
    @test navorder(t, 3) == [3]          # an empty line can still be selected
    @test isempty(orderedlines(t, 3))    # but is not written out
    @test isnothing(referenceshifts(t, 1))

    c = copy(t)
    setpick!(c, 2, Pick(2, 3.2, 1.0))
    push!(c.multiplets[2].lines, 99)
    @test npicks(t) == 1 && linesof(t, 2) == [2]

    @test relabel!(t, 2, " MO,PS ")
    @test multipletlabels(t)[2] == "MOPS"
    @test relabel!(t, 2, "Tris_HCl")
    @test multipletlabels(t)[2] == "Tris-HCl"   # an underscore would split a CSV header
    @test relabel!(t, 2, "MOPS")
    @test !relabel!(t, 2, "DSS")
    @test !relabel!(t, 2, "  ")

    erase!(q -> q.shift > 3.0, t, 2)
    @test npicks(t) == 0
end

@testset "multiplets" begin
    s = makeseries()
    t = picked(s; group=false)
    @test multipletlabels(t) == ["DSS", "MES", "MOPS", "peak3"]
    # joining keeps the current line's multiplet and its label
    @test group!(t, 3, 4) == 3
    @test multipletlabels(t) == ["DSS", "MES", "MOPS"]
    @test orderedlines(t, 3) == [3, 4]   # decreasing shift
    for i in eachindex(s)
        @test centre(t, 3, i) ≈ doublet(i) atol = 1e-5
        @test centre(t, 2, i) ≈ singlet(i) atol = 1e-3
    end
    @test referenceshifts(t, 4) ≈ dss.(1:4) atol = 1e-3

    # a plane missing one line has no centre, rather than the other line's shift
    delete!(t.picks, (4, 2))
    @test isnothing(centre(t, 3, 2))
    @test !isnothing(centre(t, 3, 1))
    # but a line never picked at all does not count
    push!(t.multiplets[3].lines, 99)
    @test !isnothing(centre(t, 3, 1))
    @test navorder(t, 3) == [3, 4, 99]

    # the reference is never absorbed into another multiplet
    @test group!(t, 2, 1) == REF
    @test multipletlabels(t) == ["DSS", "MOPS"]
    @test sort(linesof(t, REF)) == [1, 2]

    # ungrouping puts the line in a new multiplet straight after
    @test ungroup!(t, 2) == 2
    @test multipletlabels(t) == ["DSS", "peak1", "MOPS"]
    @test linesof(t, 2) == [2]
    @test ungroup!(t, 2) == 2   # already on its own

    # deleting the last line removes a multiplet, but the reference keeps an empty line
    deleteline!(t, 2)
    @test multipletlabels(t) == ["DSS", "MOPS"]
    deleteline!(t, 1)
    @test multipletlabels(t) == ["DSS", "MOPS"]
    @test linesof(t, REF) == [1]
    @test !haspicks(t, 1)
end

@testset "first-order phase outliers" begin
    @test p1outliers([0.0, 1.0])[2] == Int[]   # too few to judge
    med, ks = p1outliers([2.0, -1.0, 0.5, 1.0, 95.0, 0.0])
    @test med == 0.75
    @test ks == [5]
    # within the floor, however tight the rest of the series
    @test isempty(p1outliers([0.0, 0.0, 0.0, P1TOLERANCE - 1])[2])
    # a series that genuinely scatters widens the tolerance
    @test isempty(p1outliers([0.0, 30.0, -30.0, 60.0, -60.0])[2])
end

# ---- the brush ------------------------------------------------------------------

@testset "brush" begin
    scale = Vec2f(100, 100)   # pixels per ppm and per row
    band = densify([Point2f(0, 0), Point2f(0.1, 0)], 2, scale)
    @test length(band) == 6 && last(band) == Point2f(0.1, 0)
    @test all(diff(first.(band)) .≈ 0.02)

    # two points within a radius of each other paint one merged span
    spans = painted([Point2f(1.0, 0.5), Point2f(1.05, 0.5)], 8, scale, 0.0, 1.0)
    @test length(spans) == 1
    @test all(isapprox.(only(spans), (0.92, 1.13); atol=1e-6))
    # a point beyond the radius from the strip paints nothing
    @test isempty(painted([Point2f(1.0, 2.0)], 8, scale, 0.0, 1.0))

    s = makeseries()
    v = StackView(s)
    stroke = [Point2f(3.075, 0.1), Point2f(3.075, 3.9)]   # through the singlet in every row
    qs = track(s, v, stroke, 8, scale)
    @test [q.plane for q in qs] == 1:4
    @test all(abs(q.shift - singlet(q.plane)) < 1e-3 for q in qs)
    @test isempty(track(s, v, stroke, 8, scale; snr=200))   # every peak below the threshold
    # a stroke through only the first two rows picks only those
    @test length(track(s, v, [Point2f(3.025, 0.1), Point2f(3.025, 1.9)], 8, scale)) == 2

    # aligned offsets shift what the brush sees, not what is stored
    v.offsets[] = fill(1.0, 4)
    qs = track(s, v, [Point2f(4.075, 0.1), Point2f(4.075, 3.9)], 8, scale)
    @test all(abs(q.shift - singlet(q.plane)) < 1e-3 for q in qs)
    v.offsets[] = zeros(4)

    t = picked(s)
    q = pickat(t, 2, 3)
    @test hit(t, v, location(v, q), scale, 5) == (2, 3)
    @test hit(t, v, location(v, q), scale, 5, 1) != (2, 3)   # only line 1 considered
    @test isnothing(hit(t, v, location(v, q) .+ Point2f(1, 0), scale, 5))

    low = StackView(s; gain=0.5)
    @test nearestrow(s, low, Point2f(singlet(2), 1.2), 0.02) == 2

    others, siblings, chosen, bars = overlay(t, 3, v.gain[], v.offsets[], v.spacing, v.norm)
    @test length(chosen.points) == 4 && length(chosen.segments) == 6
    @test length(siblings.points) == 4          # the doublet's other line
    @test length(others.points) == 8            # the reference and the singlet
    @test length(others.colours) == length(others.points)
    @test length(bars.segments) == 4 * 4        # a bar and a tick in each plane
    @test length(bars.path) == 3 * 2            # joining the centres of adjacent planes
end

# ---- tables and files -------------------------------------------------------------

@testset "coordinates" begin
    @test coordinatetable(nothing, 3) == NamedTuple()
    @test coordinatetable([1, 2, 3], 3) == (; x=[1, 2, 3])
    @test coordinatetable("pH" => 6:8, 3) == (; pH=[6, 7, 8])
    @test coordinatetable((; pH=[6, 7], T=[298, 310]), 2) == (; pH=[6, 7], T=[298, 310])
    @test_throws DimensionMismatch coordinatetable([1, 2], 3)
    @test callcoordinates(nothing) === nothing
    @test callcoordinates("pH" => 6:7) == (; pH=[6, 7])
    @test exptfolder("11/pdata/1") == "11"
    @test exptfolder(joinpath("data", "11", "pdata", "2")) == joinpath("data", "11")
    @test exptfolder("11") == "11"
end

@testset "tables" begin
    s = makeseries()
    p = tracker(s, picked(s))
    header, rows = seriestable(p)
    @test header == ["source", "label", "line", "plane", "pH", "shift (ppm)", "height", "snr",
                     "centre (ppm)", "referenced (ppm)"]
    @test length(rows) == 4 + 4 + 2 * 4
    mops = filter(r -> r[2] == "MOPS", rows)
    @test [r[3] for r in mops] == [fill("1", 4); fill("2", 4)]
    @test all(parse(Float64, mops[i][6]) > parse(Float64, mops[i + 4][6]) for i in 1:4)
    # the centre is repeated down the multiplet's lines
    @test mops[2][9] == mops[6][9]
    mes = filter(r -> r[2] == "MES", rows)
    for (i, r) in enumerate(mes)
        @test parse(Float64, r[10]) ≈ singlet(i) - dss(i) atol = 2e-3
        @test r[5] == string([6.0, 6.5, 7.0, 7.5][i])
    end

    # a missing pick is NA, and so is its multiplet's centre in that plane
    delete!(tracks(p).picks, (4, 3))
    header, rows = seriestable(p)
    mops = filter(r -> r[2] == "MOPS", rows)
    @test mops[7][6] == "NA" && mops[7][9] == "NA" && mops[3][9] == "NA"

    header, rows = peaklisttable(p)
    @test header == ["label", "reference", "line", "plane", "shift (ppm)"]
    @test length(rows) == npicks(tracks(p))
    @test all(r -> (r[2] == "true") == (r[1] == "DSS"), rows)

    header, cells = widetable(p)
    @test header == ["plane", "pH", "DSS (reference)", "MES", "MOPS"]
    @test size(cells) == (4, 5)
    @test cells[3, 5] == ""   # no centre where a line is missing

    # without coordinates there is no column for them, not an empty one
    bare = tracker(s, tracks(p); coordinates=nothing)
    header, cells = widetable(bare)
    @test header == ["plane", "DSS (reference)", "MES", "MOPS"]
    @test size(cells) == (4, 4)
    @test occursin("DSS (reference)", sprint(show, MIME"text/plain"(), bare))

    # the clipboard holds the referenced shifts once there is a reference, and the measured
    # ones before then
    @test isreferenced(bare)
    text = clipboardtext(bare)
    @test startswith(text, "plane\tDSS (ppm)\tMES (ppm)\tMOPS[1] (ppm)\tMOPS[2] (ppm)\n")
    @test length(split(chomp(text), '\n')) == 5
    @test parse(Float64, split(split(text, '\n')[2], '\t')[2]) ≈ 0.0 atol = 1e-12
    unref = copy(tracks(p))
    erase!(_ -> true, unref, 1)
    @test !isreferenced(tracker(s, unref))
    @test parse(Float64, split(split(clipboardtext(tracker(s, unref)), '\n')[2], '\t')[3]) ≈
          singlet(1) atol = 1e-3

    # the wide files: one column per line, a multiplet's lines numbered from the highest shift
    header, rows = linetable(p, q -> q.shift, "ppm")
    @test header == ["plane", "pH", "DSS (ppm)", "MES (ppm)", "MOPS[1] (ppm)", "MOPS[2] (ppm)"]
    @test length(rows) == 4 && all(length.(rows) .== 6)
    @test parse(Float64, rows[2][4]) ≈ singlet(2) atol = 1e-3
    @test parse(Float64, rows[1][5]) > parse(Float64, rows[1][6])
    @test rows[3][6] == "NA"                           # the pick deleted above
    header, rows = linetable(p, q -> referencedshift(p, q), "ppm")
    @test parse(Float64, rows[2][4]) ≈ singlet(2) - dss(2) atol = 2e-3
    @test parse(Float64, rows[2][3]) ≈ 0.0 atol = 1e-12  # the reference at refshift
    header, rows = linetable(p, q -> q.height)
    @test header[end] == "MOPS[2]"                     # a height has no unit
    @test parse(Float64, rows[1][4]) ≈ 100 atol = 1
    bare = tracker(s, tracks(p); coordinates=nothing)
    @test first(linetable(bare, q -> q.shift, "ppm"))[1:2] == ["plane", "DSS (ppm)"]
    nodss = copy(tracks(p))
    erase!(_ -> true, nodss, 1)
    @test isnothing(referencedshift(tracker(s, nodss), pickat(nodss, 2, 1)))

    rows = planeresults(p)
    @test length(rows) == 16
    @test keys(first(rows)) ==
          (:label, :line, :plane, :pH, :shift, :height, :snr, :centre, :referenced)
    @test ismissing(rows[end - 1].shift)

    @test occursin("MOPS (2 lines)", sprint(show, p))
    @test occursin("DSS (reference)", sprint(show, MIME"text/plain"(), p))
end

@testset "peak list round trip" begin
    s = makeseries()
    p = tracker(s, picked(s))
    addmultiplet!(tracks(p), "empty")
    t = restore(s, peaklisttable(p)...)
    @test multipletlabels(t) == ["DSS", "MES", "MOPS", "empty"]
    @test length(orderedlines(t, 3)) == 2
    @test isempty(orderedlines(t, 4)) && length(linesof(t, 4)) == 1
    for m in 1:3, i in eachindex(s)
        @test centre(t, m, i) ≈ centre(tracks(p), m, i)
    end

    mktempdir() do dir
        path = writetable(joinpath(dir, "peaklist.csv"), ["a comment"], peaklisttable(p)...)
        header, rows = readtable(path)
        @test header == first(peaklisttable(p))
        @test length(rows) == npicks(tracks(p)) + 1
        t = restore(s, path)
        @test npicks(t) == npicks(tracks(p))

        # without a reference row, an empty reference is added in front
        header, rows = peaklisttable(p)
        rows = [[r[1], "false", r[3:end]...] for r in rows if r[1] != "DSS"]
        t = restore(s, header, rows; reference="TSP")
        @test multipletlabels(t) == ["TSP", "MES", "MOPS", "empty"]
        @test !haspicks(t, only(linesof(t, REF)))

        # picks for planes beyond the data are skipped
        short = SpectrumSeries(δaxis, s.intensities[:, 1:2], ones(2), ["a", "b"], [1, 1])
        @test_logs (:warn,) match_mode = :any restore(short, path)
    end
end

@testset "save" begin
    s = makeseries()
    p = tracker(s, picked(s))
    mktempdir() do dir
        folder = saveanalysis(p, joinpath(dir, "out"))
        for file in ["summary.txt", "peaklist.csv", "series.csv", "shifts.csv",
                     "referencedshifts.csv", "heights.csv", "stack.pdf",
                     joinpath("peaks", "MOPS.csv"), joinpath("peaks", "MOPS.pdf")]
            @test isfile(joinpath(folder, file))
        end
        @test !isfile(joinpath(folder, "processing.csv"))   # processed data were given
        @test !isfile(joinpath(folder, "results.csv"))      # nothing is fitted per peak
        summary = read(joinpath(folder, "summary.txt"), String)
        @test occursin("MOPS", summary) && occursin("2 lines", summary)
        @test !occursin("Reproduce", summary)   # no call was recorded
        @test occursin("not integrals", read(joinpath(folder, "heights.csv"), String))
        saveanalysis(p, joinpath(dir, "out"))
        @test isdir(joinpath(dir, "out_previous"))

        # without a reference there is nothing to reference to
        erase!(_ -> true, tracks(p), 1)
        folder = saveanalysis(p, joinpath(dir, "noref"))
        @test isfile(joinpath(folder, "shifts.csv"))
        @test !isfile(joinpath(folder, "referencedshifts.csv"))
    end
end
