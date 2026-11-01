# PeakTrack1D: the picking model, peak finding, brush geometry and output tables, with no
# window and no data files. Spectra are built in memory as a `SpectrumSeries`, the seam
# below which nothing touches NMRData. Processing and phasing need real FIDs and are left to
# tests alongside example data.

using NMRAnalysis
using NMRAnalysis.PeakTrack1D: SpectrumSeries, Pick, Multiplet, Tracks, REF, StackView,
                               PeakTracker
using NMRAnalysis.PeakTrack1D: climb, vertex, snap, exact, reheight!, snr, noiselevel,
                               displayscale,
                               sameaxis, planesfrom
using NMRAnalysis.PeakTrack1D: cleanlabel, multipletof, multipletlabels, linesof,
                               orderedlines, navorder,
                               centre, linecentre, referenceshifts, vacantname,
                               addmultiplet!,
                               setpick!, erase!, group!, ungroup!, deleteline!, relabel!,
                               haspicks, picksof, npicks, pickat, tracks
using NMRAnalysis.PeakTrack1D: densify, painted, track, overlay, hit, nearestrow, location
using NMRAnalysis.PeakTrack1D: coordinatetable, seriestable, peaklisttable, clipboardtext,
                               isreferenced, widetable,
                               linetable, referencedshift, summaryseries,
                               restore, readtable, saveanalysis, exptfolder,
                               callcoordinates,
                               p1outliers, P1TOLERANCE, neighbourphase, p0outliers,
                               wrap180, fixphases!, makeexternal!
using NMRAnalysis.PeakTrack1D: PickedReference, IndirectReference, referencing, alignment,
                               allreferenced, hertz, xiratio, referencename, spectralist,
                               ppmtext
using NMRTools: F19, P31, H1
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
                lorentz(80, doublet(i) + 0.02) .+ lorentz(80, doublet(i) - 0.02)
                for i in 1:n])
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

function tracker(s, t; coordinates="pH" => [6.0, 6.5, 7.0, 7.5][1:length(s)], refshift=0.0,
                 reference=PickedReference(refshift))
    fig = Figure()
    ax = Axis(fig[1, 1])
    return PeakTracker(s, coordinatetable(coordinates, length(s)), Observable(t), nothing,
                       reference, Observable(false), StackView(s), fig, ax, nothing,
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
    @test displayscale(s) == maximum(abs, s.intensities)
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
    @test vertex(x, [2.0, 2.0, 2.0], 2) == (2.0, 2.0)   # a flat top, not 0/0

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

    # heights follow the data when it is reprocessed, shifts stay where they were placed
    t = Tracks(["DSS"])
    setpick!(t, 1, snap(s, 1, dss(1), 0.01))
    q = pickat(t, 1, 1)
    scaled = SpectrumSeries(δaxis, 2 .* s.intensities, ones(4), fill("synthetic", 4), 1:4)
    reheight!(t, scaled)
    @test pickat(t, 1, 1).shift == q.shift
    @test pickat(t, 1, 1).height ≈ 2 * exact(s, 1, q.shift).height
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

    # a spectrum that couldn't be phased takes its neighbours' zero-order phase, averaged
    # round the circle
    phases = Any[(10.0, 0.0), nothing, (30.0, 0.0), nothing, (350.0, 0.0)]
    @test neighbourphase(phases, [1, 3, 5], 2) ≈ 20
    @test neighbourphase(phases, [1, 3, 5], 4) ≈ 10 atol = 1e-9
    @test neighbourphase(phases, [1, 3, 5], 6) ≈ -10

    # zero-order phases from a real ³¹P series, one of which settled ~150° from the rest
    p0 = [-155.0, -150.4, -155.8, -149.8, -156.9, -156.5, -161.8, -156.6, -157.5, -153.4,
          -161.4, -146.4, -158.3, -150.2, -163.2, -157.9, -156.6, -158.4, -156.7, -154.5,
          -160.3, -163.1, -158.2, -160.6, -160.8, 52.3, -155.2, -161.6, -163.9, -162.3,
          -161.4, -161.3]
    centre, ks = p0outliers(p0)
    @test ks == [26]
    @test -160 < centre < -155
    # deviations are measured round the circle, so phases either side of ±180° agree
    @test isempty(p0outliers([179.0, -179.0, 178.0, -178.0, 180.0])[2])
    @test p0outliers([179.0, -179.0, 178.0, -178.0, 0.0])[2] == [5]
    @test isempty(p0outliers([0.0, 90.0])[2])   # too few to judge
    @test wrap180(190.0) ≈ -170 && wrap180(-190.0) ≈ 170

    # outliers are repaired before a failed experiment borrows its neighbours' phase, so it
    # can't inherit an outlier's
    phases = Any[(0.0, 0.0), (2.0, 0.0), (-1.0, 0.0), (160.0, 0.0), nothing]
    phasing = [:auto, :auto, :auto, :auto, :failed]
    fixed, how = fixphases!(phases, phasing, fill(nothing, 5);
                            refit=(fid, p1; kw...) -> nothing)
    @test all(abs(wrap180(p[1])) < 5 for p in fixed)
    @test how == [:auto, :auto, :auto, :neighbours, :neighbours]
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

@testset "indirect referencing" begin
    s = makeseries()
    ΞF = 0.94094011
    bfH = 600.13e6
    @test xiratio(F19) == ΞF && xiratio(H1) == 1.0 && xiratio(P31) ≈ 0.40480742
    @test isnan(xiratio(nothing))
    @test hertz(600.13) == 600.13e6 && hertz(600.13e6) == 600.13e6 && isnan(hertz(nothing))

    # a 1H analysis with DSS picked exactly where it is in each plane
    hs = SpectrumSeries(δaxis, s.intensities, ones(4), fill("h", 4), 1:4; bf=bfH, xi=1.0)
    t = Tracks(["DSS"])
    foreach(i -> setpick!(t, 1, exact(hs, i, dss(i))), 1:4)
    h = tracker(hs, t)
    @test alignment(h) ≈ -dss.(1:4)
    @test referencedshift(h, 1, 3.0) ≈ 3.0 - dss(1)

    # a 19F series of the same samples whose axis is 3 ppm off the unified scale
    fs = SpectrumSeries(δaxis, s.intensities, ones(4), fill("f", 4), 1:4; nucleus="19F",
                        bf=ΞF * bfH * (1 + 3e-6), xi=ΞF)
    r = IndirectReference(h, fs)
    @test r.label == "DSS"
    @test r.frequencies ≈ [bfH * (1 + 1e-6dss(i)) for i in 1:4]
    f = tracker(fs, Tracks(["DSS"]); reference=r)
    δ′ = -60.0
    exact19F(i) = 1e6 * ((1 + 3e-6) * (1 + 1e-6δ′) / (1 + 1e-6dss(i)) - 1)
    # (ν/ν₀ − 1) × 10⁶ loses about six digits to cancellation, leaving ~1e-10 ppm
    @test referencedshift(f, 2, δ′) ≈ exact19F(2) atol = 1e-8
    @test referencedshift(f, 2, δ′) ≈ δ′ + 3 - dss(2) atol = 1e-3   # to first order
    @test allreferenced(f) && isreferenced(f)
    @test alignment(f) ≈ [referencedshift(f, i, 5.0) - 5.0 for i in 1:4]

    # a plane without DSS in the 1H analysis is not referenced in the 19F one
    t2 = copy(t)
    delete!(t2.picks, (1, 3))
    r2 = @test_logs (:warn,) match_mode = :any IndirectReference(tracker(hs, t2), fs)
    @test isnan(r2.frequencies[3])
    f2 = tracker(fs, Tracks(["DSS"]); reference=r2)
    @test isnothing(referencedshift(f2, 3, 1.0)) && !allreferenced(f2)
    @test isnothing(alignment(f2))

    # spectra that don't match one to one, or differ in field, are refused
    three = SpectrumSeries(δaxis, s.intensities[:, 1:3], ones(3), fill("f", 3), 1:3;
                           bf=ΞF * bfH, xi=ΞF)
    @test_throws DimensionMismatch IndirectReference(h, three)
    @test_throws DimensionMismatch SpectrumSeries(δaxis, s.intensities, ones(4),
                                                  fill("f", 4),
                                                  1:4; bf=[bfH], xi=ΞF)

    # the field is judged from BF1, so an axis re-referenced by 300 ppm still matches
    h1 = tracker(SpectrumSeries(δaxis, s.intensities, ones(4), fill("h", 4), 1:4; bf=bfH,
                                bf1=bfH, xi=1.0), t)
    shifted = SpectrumSeries(δaxis, s.intensities, ones(4), fill("f", 4), 1:4;
                             bf=ΞF * bfH * (1 + 300e-6), bf1=ΞF * bfH, xi=ΞF)
    @test IndirectReference(h1, shifted) isa IndirectReference
    elsewhere = SpectrumSeries(δaxis, s.intensities, ones(4), fill("f", 4), 1:4;
                               bf=ΞF * bfH, bf1=ΞF * 500.13e6, xi=ΞF)
    @test_throws ArgumentError IndirectReference(h1, elsewhere)
    otherfield = SpectrumSeries(δaxis, s.intensities, ones(4), fill("f", 4), 1:4;
                                bf=ΞF * 500.13e6, xi=ΞF)
    @test_throws ArgumentError IndirectReference(h, otherfield)
    @test_throws ArgumentError IndirectReference(f, fs)                         # not picked
    @test_throws ArgumentError IndirectReference(tracker(fs, Tracks(["X"])), fs) # not 1H
    @test_throws ArgumentError IndirectReference(h, makeseries())               # no Ξ

    # a reference picked in an earlier peak list stays a peak of its own when the reference
    # comes from elsewhere, rather than having its picks relabelled
    local3 = Tracks(["TFA"])
    setpick!(local3, 1, Pick(1, -76.5, 1.0))
    addmultiplet!(local3, "peak1")
    makeexternal!(local3, "DSS")
    @test multipletlabels(local3) == ["DSS", "TFA", "peak1"]
    @test !haspicks(local3, only(linesof(local3, REF)))
    @test pickat(local3, 1, 1).shift == -76.5 && multipletof(local3, 1) == 2
    empty3 = makeexternal!(Tracks(["reference"]), "DSS")
    @test multipletlabels(empty3) == ["DSS"]

    # the reference keyword
    @test referencing("D_SS", 0.5, fs) == ("D-SS", PickedReference(0.5))
    @test first(referencing(h, 0.0, fs)) == "DSS"
    @test referencing(("DSS", r), 0.0, fs) == ("DSS", r)

    # saved output says how it was referenced
    t3 = Tracks(["DSS"])
    setpick!(t3, addmultiplet!(t3, "TFA"), exact(fs, 1, δ′))
    mktempdir() do dir
        folder = saveanalysis(tracker(fs, t3; reference=r), joinpath(dir, "out"))
        text = read(joinpath(folder, "referencedshifts.csv"), String)
        @test occursin("indirectly", text) && occursin("TFA (ppm)", text)
        @test occursin("indirectly", read(joinpath(folder, "summary.txt"), String))
    end
end

@testset "message wording" begin
    @test referencename(Tracks(["reference"])) == "the reference peak"
    @test referencename(Tracks(["DSS"])) == "DSS"
    @test spectralist([3], 32) == "spectrum 3"
    @test spectralist([3, 7, 12], 32) == "spectra 3, 7 and 12"
    @test spectralist(1:20, 32) == "20 of the 32 spectra"
    @test ppmtext(0.0) == "0 ppm" && ppmtext(4.7) == "4.7 ppm"
end

@testset "summary panel" begin
    s = makeseries()
    p = tracker(s, picked(s))
    d = summaryseries(p, 3)                        # a line of the doublet
    @test d.m == 3 && d.label == "MOPS"
    @test d.xlabel == "pH"
    @test startswith(d.ylabel, "shift referenced to DSS")
    @test length(d.centre) == 4
    @test d.centre[2][2] ≈ doublet(2) - dss(2) atol = 2e-3
    @test length(d.lines) == 2 * (4 + 1)           # two lines, each ending in a gap
    @test isnan(d.lines[5][1])
    @test d.heights[1][2] ≈ 80 atol = 6

    # a missing pick leaves a gap rather than joining its neighbours
    delete!(tracks(p).picks, (4, 3))
    d = summaryseries(p, 3)
    @test isnan(d.centre[3][2]) && isnan(d.lines[8][2])

    # without the reference in every plane, shifts are as measured
    delete!(tracks(p).picks, (1, 2))
    d = summaryseries(p, 2)
    @test startswith(d.ylabel, "shift as measured")
    @test d.centre[1][2] ≈ singlet(1) atol = 1e-3

    # the vacant slot shows nothing, with the same fields
    e = summaryseries(p, 0)
    @test e.m == 0 && isempty(e.centre) && keys(e) == keys(d)

    # points run in order of the coordinate, not of the planes
    q = tracker(s, tracks(p); coordinates="T" => [310, 300, 290, 280])
    @test issorted(first.(summaryseries(q, 2).centre))
    @test summaryseries(tracker(s, tracks(p); coordinates=nothing), 2).xlabel == "spectrum"
end

@testset "tables" begin
    s = makeseries()
    p = tracker(s, picked(s))
    header, rows = seriestable(p)
    @test header ==
          ["source", "label", "line", "plane", "pH", "shift (ppm)", "height", "snr",
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
    @test header == ["spectrum", "pH", "DSS (reference)", "MES", "MOPS"]
    @test size(cells) == (4, 5)
    @test cells[3, 5] == ""   # no centre where a line is missing

    # without coordinates there is no column for them, not an empty one
    bare = tracker(s, tracks(p); coordinates=nothing)
    header, cells = widetable(bare)
    @test header == ["spectrum", "DSS (reference)", "MES", "MOPS"]
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
    @test header ==
          ["plane", "pH", "DSS (ppm)", "MES (ppm)", "MOPS[1] (ppm)", "MOPS[2] (ppm)"]
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
        touch(joinpath(folder, "figure-20261101-120000.pdf"))
        saveanalysis(p, joinpath(dir, "out"))
        @test isdir(joinpath(dir, "out_previous"))
        # figures saved with f are copied into the new folder, and kept in the old one
        @test isfile(joinpath(dir, "out", "figure-20261101-120000.pdf"))
        @test isfile(joinpath(dir, "out_previous", "figure-20261101-120000.pdf"))

        # without a reference there is nothing to reference to
        erase!(_ -> true, tracks(p), 1)
        folder = saveanalysis(p, joinpath(dir, "noref"))
        @test isfile(joinpath(folder, "shifts.csv"))
        @test !isfile(joinpath(folder, "referencedshifts.csv"))
    end
end
