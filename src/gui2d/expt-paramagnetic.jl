# Paramagnetic relaxation enhancement, shared by pre2d and lpre2d.
#
# These experiments do not use the shared fit of fitting.jl, which gives every plane a free
# amplitude: the PRE is measured precisely by how the amplitude changes from plane to plane.
# Each peak instead has one amplitude, its diamagnetic height, and a model giving its
# linewidths and relative height in every plane from a few parameters (`peakmodel`). The
# amplitudes are still solved for linearly at every step (variable projection, as in
# fitting.jl), and the model's own parameters join the positions and linewidths in the
# nonlinear fit.
#
# A concrete experiment provides `modelparameters`, `planecondition` and `peakmodel`, and
# has the fields `concentrations`, `paired`, `expttype`, `Trelax` and `skipplanes`.

abstract type ParamagneticExperiment <: FixedPeakExperiment end

"""
    modelparameters(expt) -> Vector{NamedTuple}

The parameters fitted for every peak besides its position and diamagnetic linewidths, each
as `(key, label, initial, lower, upper, scale)`. `scale` is the size of a typical value, so
that a rate per unit concentration is fitted on the same footing whatever the unit.
"""
function modelparameters end

# ---- the physics ---------------------------------------------------------------

# One-bond couplings setting the nominal transfer delays: amide ¹H-¹⁵N and methyl ¹H-¹³C
const ONEBONDJ = Dict(N15 => 92.0, C13 => 125.0)

"""
    transfertime(Trelax, spec) -> Float64

`Trelax` in seconds, or for `:auto` the 1/¹J for the heteronucleus in F2 of `spec`: two
INEPTs of 1/(2J) in an HSQC, or the two 1/(2J) delays either side of t₁ in an HMQC.
"""
transfertime(Trelax::Real, spec) = Float64(Trelax)
function transfertime(Trelax::Symbol, spec)
    Trelax === :auto ||
        throw(ArgumentError("Trelax must be a time in seconds or :auto, " *
                            "not $(repr(Trelax))"))
    nuc = spec[F2Dim, :nucleus]
    haskey(ONEBONDJ, nuc) ||
        throw(ArgumentError("cannot set Trelax automatically for a $(label(spec, F2Dim)) " *
                            "indirect dimension: give it in seconds"))
    return 1 / ONEBONDJ[nuc]
end

"""
    transverseplane(expt, c, R2x, R2y, Γ2) -> (R2x, R2y, height)

Linewidths and relative peak height at paramagnetic concentration `c`. Γ₂c adds to the ¹H
linewidth (and in an HMQC to the multiple-quantum linewidth in F2), and the height falls by
the broadening and by exp(-Γ₂c Trelax) during the transfer delays, giving
I_para/I_dia = R₂ exp(-Γ₂t)/(R₂ + Γ₂) for an HSQC (Battiste & Wagner, Biochemistry 39,
5355 (2000); Iwahara, Tang & Clore, J. Magn. Reson. 184, 185 (2007)).
"""
function transverseplane(expt::ParamagneticExperiment, c, R2x, R2y, Γ2)
    Γ = Γ2 * c
    R2xp = R2x + Γ
    R2yp = expt.expttype == :hmqc ? R2y + Γ : R2y
    return R2xp, R2yp, exp(-Γ * expt.Trelax) * (R2x / R2xp) * (R2y / R2yp)
end

# The fitted parameters of a peak, in the order the fit holds them
function fitparamnames(expt::ParamagneticExperiment)
    return [:x, :y, :R2x0, :R2y0,
            (m.key for m in modelparameters(expt))...]
end
npeakparams(expt::ParamagneticExperiment) = 4 + length(modelparameters(expt))
fitvalues(peak, expt) = [peak.parameters[s].value[][1] for s in fitparamnames(expt)]

"Linewidths and relative height in plane `i` of a peak with fit parameters `φp`."
function planeshape(expt::ParamagneticExperiment, i, φp)
    return peakmodel(expt, planecondition(expt, i), φp[3], φp[4], φp[5:end])
end

# ---- planes and references -----------------------------------------------------

"Check the sequence type and that the concentrations give one per spectrum, and two or more."
function checksamples(expttype, concentrations, n, skipplanes)
    expttype in (:hsqc, :hmqc) ||
        throw(ArgumentError("type must be :hsqc or :hmqc, not $(repr(expttype))"))
    length(concentrations) == n ||
        throw(ArgumentError("got $(length(concentrations)) concentrations for $n spectra"))
    all(≥(0), concentrations) ||
        throw(ArgumentError("concentrations must not be negative"))
    used = [c for (i, c) in enumerate(concentrations) if i ∉ skipplanes]
    length(unique(used)) ≥ 2 ||
        throw(ArgumentError("PRE needs spectra at two or more concentrations"))
    return nothing
end

usedplanes(expt) = [i for i in 1:nslices(expt) if i ∉ skipset(expt)]
paramagneticplanes(expt) = [i for i in usedplanes(expt) if expt.concentrations[i] > 0]

"""
    relativeintensity(peak, expt, i) -> Measurement or nothing

The height of `peak` in paramagnetic plane `i` relative to its mean height in the
diamagnetic planes it is compared with (see `referenceplanes`), or `nothing` where there are
none.
"""
function relativeintensity(peak, expt::ParamagneticExperiment, i)
    refs = referenceplanes(expt, i)
    isempty(refs) && return nothing
    a = peak.parameters[:amp].value[] .± peak.parameters[:amp].uncertainty[]
    return a[i] / (sum(a[refs]) / length(refs))
end

# A single paramagnetic plane gives one ratio per peak, worth a column in results.csv
hasratio(expt::ParamagneticExperiment) = length(paramagneticplanes(expt)) == 1

# ---- data and peaks ------------------------------------------------------------

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, zlabels, ::Type{<:ParamagneticExperiment})
    @debug "Preparing spec data for PRE experiment: $inputfilenames"
    spectra = map(loadnmr, inputfilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)
    return SpecData(spectra, x, y, z ./ σ[1], σ ./ σ[1], collect(String, zlabels))
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::ParamagneticExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    peak = Peak(initialposition, label, xradius, yradius)
    n = nslices(expt)
    x0, y0 = initialposition
    amp0 = map(1:n) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end

    # what the fit varies: diamagnetic linewidths, the model, one amplitude
    peak.parameters[:R2x0] = Parameter("R2x0", 20.0; minvalue=1.0, maxvalue=100.0)
    peak.parameters[:R2y0] = Parameter("R2y0", 15.0; minvalue=1.0, maxvalue=100.0)
    for m in modelparameters(expt)
        peak.parameters[m.key] = Parameter(m.label, m.initial; minvalue=m.lower,
                                           maxvalue=m.upper)
    end
    peak.parameters[:amp0] = Parameter("Diamagnetic amplitude", maximum(amp0))
    # what follows in each plane: linewidths, measured and fitted heights
    peak.parameters[:R2x] = Parameter("R2x", fill(20.0, n))
    peak.parameters[:R2y] = Parameter("R2y", fill(15.0, n))
    peak.parameters[:amp] = Parameter("Amplitude", amp0)
    peak.parameters[:ampfit] = Parameter("Fitted amplitude", copy(amp0))

    for m in modelparameters(expt)
        peak.postparameters[m.key] = Parameter(m.label, 0.0)
    end
    peak.postparameters[:R2x0] = Parameter("R2x0", 0.0)
    peak.postparameters[:R2y0] = Parameter("R2y0", 0.0)
    hasratio(expt) && (peak.postparameters[:Irel] = Parameter("I/Iref", 0.0))

    push!(expt.peaks[], peak)
    return notify(expt.peaks)
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::ParamagneticExperiment)
    for key in keys(peak.postparameters)
        key == :Irel && continue
        peak.postparameters[key].value[] .= peak.parameters[key].value[]
        peak.postparameters[key].uncertainty[] .= peak.parameters[key].uncertainty[]
    end
    r = hasratio(expt) ? relativeintensity(peak, expt, only(paramagneticplanes(expt))) :
        nothing
    if !isnothing(r)
        peak.postparameters[:Irel].value[] .= Measurements.value(r)
        peak.postparameters[:Irel].uncertainty[] .= Measurements.uncertainty(r)
    end
    return peak.postfitted[] = true
end

# the fit, not the measured heights, so the displayed fit is the model
function simulate!(z, peak::Peak, expt::ParamagneticExperiment)
    amp = peak.parameters[:ampfit].value[]
    for i in 1:nslices(expt)
        xi, yi, block = peakshape(PlaneGrid(expt, i), shapeparams(peak, i)...,
                                  peak.xradius[], peak.yradius[])
        z[i][xi, yi] .+= amp[i] .* block
    end
    return z
end

function fittedamplitudes(peak, expt::ParamagneticExperiment)
    peak.postfitted[] || return fill(NaN, nslices(expt))
    return collect(Float64, peak.parameters[:ampfit].value[])
end

# ---- the fit ---------------------------------------------------------------------

function fitcluster!(peaks, expt::ParamagneticExperiment, check; warm=false)
    masks = mask(peaks, expt)
    groups = fitgroups(expt, masks, 1:nslices(expt); skip=skipset(expt))
    result = fitmodel(peaks, groups, expt, check; warm)
    commitmodel!(peaks, expt, result)
    return result.status
end

"Fit coordinates for `peaks`: positions and linewidths as `ShapeFit`, then the model."
function modelshapefit(peaks, expt::ParamagneticExperiment, grid)
    origin, scale, lower, upper = Float64[], Float64[], Float64[], Float64[]
    for peak in peaks
        x0 = peak.parameters[:x].initialvalue[][1]
        y0 = peak.parameters[:y].initialvalue[][1]
        R2x, R2y = peak.parameters[:R2x0], peak.parameters[:R2y0]
        append!(origin, (x0, y0, 0.0, 0.0))
        append!(scale, (peak.xradius[], peak.yradius[], 1.0, 1.0))
        append!(lower, (-1.0, -1.0, R2x.minvalue[], R2y.minvalue[]))
        append!(upper,
                (1.0, 1.0, max(R2x.maxvalue[], maxR2(grid.xaxis, x0, peak.xradius[])),
                 max(R2y.maxvalue[], maxR2(grid.yaxis, y0, peak.yradius[]))))
        for m in modelparameters(expt)
            push!(origin, 0.0)
            push!(scale, m.scale)
            push!(lower, m.lower / m.scale)
            push!(upper, m.upper / m.scale)
        end
    end
    return ShapeFit(origin, scale, lower, upper)
end

"""
    planedesign(expt, group, i, φ, radii; scaled=true) -> Matrix

The design matrix of plane `i` of `group`: column `p` is peak `p` with the linewidths it has
in that plane, at its diamagnetic height (`scaled`) or at unit height.
"""
function planedesign(expt::ParamagneticExperiment, g, i, φ, radii; scaled=true)
    k = npeakparams(expt)
    A = zeros(count(g.mask), length(radii))
    block = zeros(length(g.grid.x), length(g.grid.y))
    for p in eachindex(radii)
        φp = φ[(k * (p - 1) + 1):(k * p)]
        R2x, R2y, s = planeshape(expt, i, φp)
        xi, yi, b = peakshape(g.grid, φp[1], φp[2], R2x, R2y, radii[p]...)
        fill!(block, 0.0)
        block[xi, yi] .= (scaled ? s : 1.0) .* b
        A[:, p] .= view(vec(block), g.mask)
    end
    return A
end

# the planes of `fitted` ((group, positions in group) pairs) stacked in turn
function modeldesign(expt, fitted, φ, radii)
    return reduce(vcat,
                  [planedesign(expt, g, g.planes[k], φ, radii)
                   for (g, ks) in fitted for k in ks])
end
observed(fitted) = reduce(vcat, [vec(g.z[:, ks]) for (g, ks) in fitted])

"""
    fitmodel(peaks, groups, expt, check; warm=false) -> NamedTuple

Fit `peaks` across every unskipped plane of `groups` with one amplitude each, returning the
fitted parameters `φ` and their errors, the amplitudes `a`, and for every peak (row) and
plane (column) its linewidths, its fitted height and its height measured at those
linewidths, each with an error. Errors come from the covariance of the full problem,
amplitudes included, σ²(JᵀJ)⁻¹.
"""
function fitmodel(peaks, groups, expt::ParamagneticExperiment, check; warm=false)
    k = npeakparams(expt)
    f = modelshapefit(peaks, expt, first(groups).grid)
    radii = [(p.xradius[], p.yradius[]) for p in peaks]
    names = fitparamnames(expt)
    φ0 = reduce(vcat,
                [[(warm ? p.parameters[s].value[] : p.parameters[s].initialvalue[])[1]
                  for s in names] for p in peaks])
    θ0 = clamp.((φ0 .- f.origin) ./ f.scale, f.lower, f.upper)
    fitted = [(g, findall(g.active)) for g in groups if any(g.active)]
    z = observed(fitted)

    function resid(θ)
        check()
        A = modeldesign(expt, fitted, physical(f, θ), radii)
        return z .- A * amplitudes(A, z)
    end

    sol = LsqFit.lmfit(resid, θ0, Float64[]; lower=f.lower, upper=f.upper,
                       autodiff=:finite, maxIter=200, x_tol=1e-8, g_tol=1e-10)
    θ = coef(sol)
    check()
    φ = physical(f, θ)
    A = modeldesign(expt, fitted, φ, radii)
    a = amplitudes(A, z)

    # Jacobian of the model A(θ)a in θ and a, for the covariance of both
    nθ = length(θ)
    J = zeros(length(z), nθ + length(a))
    for j in 1:nθ
        check()
        h = 1e-6 * max(abs(θ[j]), 1.0)
        θp, θm = copy(θ), copy(θ)
        θp[j] += h
        θm[j] -= h
        J[:, j] .= (modeldesign(expt, fitted, physical(f, θp), radii) * a .-
                    modeldesign(expt, fitted, physical(f, θm), radii) * a) ./ (2h)
    end
    J[:, (nθ + 1):end] .= A
    mse = sum(abs2, z .- A * a) / max(length(z) - size(J, 2), 1)
    C = mse .* covinverse(J' * J)
    σ = sqrt.(max.(diag(C), 0.0))

    n = nslices(expt)
    np = length(peaks)
    R2x, σR2x, R2y, σR2y = zeros(np, n), zeros(np, n), zeros(np, n), zeros(np, n)
    fit, σfit = zeros(np, n), zeros(np, n)
    for p in 1:np, i in 1:n
        idx = (k * (p - 1) + 1):(k * p)
        cols = [collect(idx); nθ + p]
        # linewidths and fitted height in plane i, and their gradient in (θp, ap)
        q = function (t)
            wx, wy, s = planeshape(expt, i, f.origin[idx] .+ f.scale[idx] .* t)
            return [wx, wy, a[p] * s]
        end
        v = q(θ[idx])
        G = zeros(3, k + 1)
        for j in 1:k
            h = 1e-6 * max(abs(θ[idx[j]]), 1.0)
            tp, tm = θ[idx], θ[idx]
            tp[j] += h
            tm[j] -= h
            G[:, j] .= (q(tp) .- q(tm)) ./ (2h)
        end
        G[3, k + 1] = planeshape(expt, i, φ[idx])[3]
        e = sqrt.(max.(diag(G * C[cols, cols] * G'), 0.0))
        R2x[p, i], R2y[p, i], fit[p, i] = v
        σR2x[p, i], σR2y[p, i], σfit[p, i] = e
    end

    # each plane's own heights, measured at the linewidths the model gives it
    amp, σamp = zeros(np, n), zeros(np, n)
    for g in groups, (j, i) in enumerate(g.planes)
        Ai = planedesign(expt, g, i, φ, radii; scaled=false)
        amp[:, i] .= amplitudes(Ai, g.z[:, j])
        σamp[:, i] .= sqrt.(max.(mse .* diag(covinverse(Ai' * Ai)), 0.0))
    end

    shape = [j for j in eachindex(θ) if mod1(j, k) ≤ 4]
    atbound = any(j -> θ[j] ≤ f.lower[j] + 1e-6 * abs(f.lower[j]) + 1e-9 ||
                       θ[j] ≥ f.upper[j] - 1e-6 * abs(f.upper[j]) - 1e-9, shape)
    status = !sol.converged ? :maxiter : atbound ? :bound : :converged
    return (; φ, σφ=f.scale .* σ[1:nθ], a, σa=σ[(nθ + 1):end], R2x, σR2x, R2y, σR2y,
            fit, σfit, amp, σamp, status)
end

"Write a model fit into `peaks`."
function commitmodel!(peaks, expt::ParamagneticExperiment, r)
    k = npeakparams(expt)
    names = fitparamnames(expt)
    setvalue!(par, i, v, e) = (par.value[][i]=v; par.uncertainty[][i]=e)
    for (p, peak) in enumerate(peaks)
        for (j, sym) in enumerate(names)
            setvalue!(peak.parameters[sym], 1, r.φ[k * (p - 1) + j], r.σφ[k * (p - 1) + j])
        end
        setvalue!(peak.parameters[:amp0], 1, r.a[p], r.σa[p])
        for i in 1:nslices(expt)
            setvalue!(peak.parameters[:R2x], i, r.R2x[p, i], r.σR2x[p, i])
            setvalue!(peak.parameters[:R2y], i, r.R2y[p, i], r.σR2y[p, i])
            setvalue!(peak.parameters[:amp], i, r.amp[p, i], r.σamp[p, i])
            setvalue!(peak.parameters[:ampfit], i, r.fit[p, i], r.σfit[p, i])
        end
    end
    return peaks
end

# ---- text and units --------------------------------------------------------------

israte(expt::ParamagneticExperiment, key) = key in (m.key for m in modelparameters(expt))

# A rate per unit concentration is per a unit that was never given
function displayunit(expt::ParamagneticExperiment, key)
    key in (:R1, :R2x0, :R2y0) && return "s⁻¹"
    return expt.paired ? "s⁻¹" : "s⁻¹ per unit concentration"
end

function paramunit(e::ParamagneticExperiment, name::Symbol)
    name in (:R1, :R2x0, :R2y0) && return "s-1"
    israte(e, name) && return e.paired ? "s-1" : ""
    return paramunit(name)
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::ParamagneticExperiment, idx)
    return "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::ParamagneticExperiment, idx)
    idx == 0 && return "No peak selected"
    peak = expt.peaks[][idx]
    peak.postfitted[] || return "Peak: $(peak.label[])\nNot fitted"
    lines = ["Peak: $(peak.label[])"]
    for m in modelparameters(expt)
        push!(lines,
              "$(m.label): $(valueerror(peak.parameters[m.key])) $(displayunit(expt, m.key))")
    end
    hasratio(expt) && push!(lines, "I/Iref: $(valueerror(peak.postparameters[:Irel]))")
    append!(lines,
            ["",
             "δX: $(valueerror(peak.parameters[:x])) ppm",
             "δY: $(valueerror(peak.parameters[:y])) ppm",
             "Diamagnetic amplitude: $(valueerror(peak.parameters[:amp0]))",
             "Diamagnetic X linewidth: $(valueerror(peak.parameters[:R2x0])) s⁻¹",
             "Diamagnetic Y linewidth: $(valueerror(peak.parameters[:R2y0])) s⁻¹"])
    return join(lines, "\n")
end

# ---- peak plots --------------------------------------------------------------------

struct ParamagneticVisualisation <: VisualisationStrategy end
visualisationtype(::ParamagneticExperiment) = ParamagneticVisualisation()

const SERIESCOLOURS = [:steelblue, :orangered, :seagreen, :purple, :goldenrod]
function seriescolour(j, m=length(SERIESCOLOURS))
    return m ≤ length(SERIESCOLOURS) ? SERIESCOLOURS[j] :
           cgrad(:viridis, m; categorical=true)[j]
end

"The planes drawn as one series in a peak plot, by concentration, as `label => planes`."
function seriesplanes(expt::ParamagneticExperiment)
    cs = sort(unique(expt.concentrations))
    return map(cs) do c
        label = c == 0 ? "Diamagnetic" : expt.paired ? "Paramagnetic" : "c = $c"
        return label => findall(==(c), expt.concentrations)
    end
end

function serieslegend!(ax, expt)
    series = seriesplanes(expt)
    m = length(series)
    elements = [MarkerElement(; color=seriescolour(j, m), marker=:circle)
                for j in 1:m]
    return axislegend(ax, elements, first.(series); framevisible=false,
                      position=legendposition(expt))
end
legendposition(::ParamagneticExperiment) = :rt

# Two concentrations are best seen as cross-sections, more as a curve through the heights
crosssectionview(expt::ParamagneticExperiment) = length(unique(expt.concentrations)) == 2
haslegend(expt::ParamagneticExperiment) = crosssectionview(expt)

function peakaxes(panel, expt::ParamagneticExperiment)
    crosssectionview(expt) || return (Axis(panel[1, 1]; curvelabels(expt)...),)
    return (Axis(panel[1, 1]; xlabel="δX / ppm", xreversed=true),
            Axis(panel[1, 2]; xlabel="δY / ppm", xreversed=true))
end

function drawpeak!(axes, peak, expt::ParamagneticExperiment)
    foreach(ax -> hlines!(ax, [0]; color=:grey50, linewidth=1), axes)
    isnothing(peak) && return nothing
    return crosssectionview(expt) ? drawcrosssections!(axes..., peak, expt) :
           drawcurve!(only(axes), peak, expt)
end

function drawcrosssections!(axX, axY, peak, expt::ParamagneticExperiment)
    series = seriesplanes(expt)
    skip = skipset(expt)
    for (j, (_, planes)) in enumerate(series), i in planes
        colour = i in skip ? :gray60 : seriescolour(j, length(series))
        x, y = expt.specdata.x[i], expt.specdata.y[i]
        x0, y0 = peak.parameters[:x].value[][i], peak.parameters[:y].value[][i]
        ix = abs.(x .- x0) .≤ peak.xradius[]
        iy = abs.(y .- y0) .≤ peak.yradius[]
        ix0, iy0 = findnearest(x, x0), findnearest(y, y0)
        z, zfit = expt.specdata.z[i], expt.specdata.zfit[][i]
        scatter!(axX, x[ix], z[ix, iy0]; color=colour, markersize=5)
        lines!(axX, x[ix], zfit[ix, iy0]; color=colour)
        scatter!(axY, y[iy], z[ix0, iy]; color=colour, markersize=5)
        lines!(axY, y[iy], zfit[ix0, iy]; color=colour)
    end
end

function plot_peak!(panel, peak, expt, ::ParamagneticVisualisation)
    axes = peakaxes(panel, expt)
    drawpeak!(axes, peak, expt)
    haslegend(expt) && !isnothing(peak) && serieslegend!(first(axes), expt)
    return axes
end

completestate!(state, expt, ::ParamagneticVisualisation) = state

function makepeakplot!(gui, state, expt, vis::ParamagneticVisualisation)
    axes = peakaxes(gui[:panelpeakplot], expt)
    if length(axes) == 2
        gui[:axpeakplotX], gui[:axpeakplotY] = axes
    else
        gui[:axpeakplot] = only(axes)
    end
    haslegend(expt) && serieslegend!(first(axes), expt)
    redraw = function (_)
        foreach(empty!, axes)
        drawpeak!(axes, state[:current_peak][], expt)
        return foreach(autolimits!, axes)
    end
    on(redraw, state[:current_peak])
    on(redraw, expt.specdata.zfit)
    return redraw(nothing)
end

# ---- summary plots -----------------------------------------------------------------

"""
    summaryplot(expt::ParamagneticExperiment; title, size, include_unassigned)

Stacked panels against residue number: the fitted rates, each peak's height in every
paramagnetic plane relative to the diamagnetic reference, and the diamagnetic relaxation
rates. Pass `param` to plot a single parameter instead, as for any other experiment.
"""
function summaryplot(expt::ParamagneticExperiment; param=nothing, title="", size=nothing,
                     include_unassigned=false, kwargs...)
    isnothing(param) ||
        return invoke(summaryplot, Tuple{Any}, expt; param, title, size,
                      include_unassigned, kwargs...)
    peaks = filter(p -> p.postfitted[], sortedpeaks(expt))
    if !include_unassigned
        assigned = filter(p -> parse_label(p.label[]).resnum > 0, peaks)
        isempty(assigned) || (peaks = assigned)
    end
    atoms = has_atom_labels(p.label[] for p in peaks)
    if atoms
        peaks = sort(peaks; by=p -> (l=parse_label(p.label[]); (l.resnum, l.atom)))
    end
    xs = atoms ? collect(1.0:length(peaks)) :
         Float64[parse_label(p.label[]).resnum for p in peaks]

    panels = summarypanels(expt, peaks)
    n = length(panels)
    fig = Figure(; size=something(size, (800, 180 * n + 80)))
    axes = Axis[]
    for (row, panel) in enumerate(panels)
        ax = Axis(fig[row, 1]; ylabel=panel.ylabel, title=(row == 1 ? title : ""),
                  xlabel=(row == n && !atoms ? "Residue number" : ""),
                  xgridvisible=false, ygridvisible=false)
        if atoms
            ax.xticks = (xs, [p.label[] for p in peaks])
            ax.xticklabelrotation = π / 4
            ax.xticklabelsvisible = row == n
        end
        drawpanel!(ax, xs, panel)
        length(panel.series) > 1 &&
            Legend(fig[row, 2], ax; framevisible=false, labelsize=11)
        push!(axes, ax)
    end
    n > 1 && linkxaxes!(axes...)
    return fig
end

"A summary panel of one parameter of `peaks`, from their fitted values."
function parampanel(peaks, key, ylabel)
    v = [p.postparameters[key].value[][1] for p in peaks]
    e = [p.postparameters[key].uncertainty[][1] for p in peaks]
    return (; ylabel, series=[(label="", values=v, errors=e)], reference=0.0)
end

function ratepanel(expt::ParamagneticExperiment, peaks, key, symbol)
    unit = expt.paired ? "s⁻¹" : "s⁻¹ per unit conc."
    return parampanel(peaks, key, "$symbol / $unit")
end

"Heights relative to the diamagnetic reference, a series for each paramagnetic plane."
function intensitypanel(expt::ParamagneticExperiment, peaks)
    series = map(paramagneticplanes(expt)) do i
        r = [something(relativeintensity(p, expt, i), NaN ± NaN) for p in peaks]
        return (label=planelabel(expt, i), values=Measurements.value.(r),
                errors=Measurements.uncertainty.(r))
    end
    return (; ylabel="I / I_ref", series, reference=1.0)
end

"The diamagnetic linewidths, a series for each dimension."
function transversepanel(expt::ParamagneticExperiment, peaks)
    spec = first(expt.specdata.nmrdata)
    nx, ny = label(spec, F1Dim), label(spec, F2Dim)
    ylab = expt.expttype == :hmqc ? "$nx-$ny MQ" : ny
    series = [(label=nx, values=[p.postparameters[:R2x0].value[][1] for p in peaks],
               errors=[p.postparameters[:R2x0].uncertainty[][1] for p in peaks]),
              (label=ylab, values=[p.postparameters[:R2y0].value[][1] for p in peaks],
               errors=[p.postparameters[:R2y0].uncertainty[][1] for p in peaks])]
    return (; ylabel="Diamagnetic R₂ / s⁻¹", series, reference=0.0)
end

function drawpanel!(ax, xs, panel)
    hlines!(ax, [panel.reference]; color=:grey50, linewidth=1)
    m = length(panel.series)
    for (j, s) in enumerate(panel.series)
        x = xs .+ (m > 1 ? 0.6 * ((j - 0.5) / m - 0.5) : 0.0)
        ok = isfinite.(s.values)
        colour = seriescolour(j, m)
        errs = [isfinite(e) ? e : 0.0 for e in s.errors[ok]]
        errorbars!(ax, x[ok], s.values[ok], errs; whiskerwidth=4, color=colour)
        scatter!(ax, x[ok], s.values[ok]; color=colour, markersize=7,
                 label=(isempty(s.label) ? nothing : s.label))
    end
end
