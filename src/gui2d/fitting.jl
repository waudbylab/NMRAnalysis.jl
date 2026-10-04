# Lineshape fitting, shared by every 2D experiment.
#
# A peak is the product of two one-dimensional lineshapes, truncated at its fitting radius
# and scaled by π²·R2x·R2y so that its amplitude is the peak height. The model is therefore
# linear in the amplitudes, which are solved for exactly at every step of the nonlinear fit
# of positions and linewidths (variable projection: Golub & Pereyra, SIAM J. Numer. Anal.
# 10, 413 (1973); O'Leary & Rust, Comput. Optim. Appl. 54, 579 (2013)). For a fixed-peak
# series this leaves four nonlinear parameters per peak however many planes there are.

"Seconds a cluster may take before its fit is stopped, keeping the window responsive."
const FIT_TIME_BUDGET = 30.0

"Seconds allowed when continuing an unfinished fit at the user's request."
const CONTINUE_TIME_BUDGET = 300.0

"Thrown inside a fit whose inputs have changed, or which the user cancelled."
struct FitCancelled <: Exception end

"Thrown inside a fit that has run past its time budget."
struct FitTimeout <: Exception end

# ---- grids and peak shapes --------------------------------------------------------

"""
    PlaneGrid

The chemical shifts of one plane within a fitting window, with the axes that convert
them to frequencies and carry the apodisation.
"""
struct PlaneGrid
    x::Vector{Float64}
    y::Vector{Float64}
    ωx::Vector{Float64}
    ωy::Vector{Float64}
    xaxis::Any
    yaxis::Any
end

function PlaneGrid(expt::Experiment, i, xrange=Colon(), yrange=Colon())
    xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
    yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
    x = collect(Float64, expt.specdata.x[i][xrange])
    y = collect(Float64, expt.specdata.y[i][yrange])
    return PlaneGrid(x, y, 2π .* hz(x, xaxis), 2π .* hz(y, yaxis), xaxis, yaxis)
end

"""
    peakshape(grid, x0, y0, R2x, R2y, xradius, yradius) -> (xi, yi, block)

A unit-amplitude peak on `grid`: `block` holds its values on the rows `xi` and columns
`yi` lying within the radii of `(x0, y0)`, and it is zero elsewhere.
"""
function peakshape(g::PlaneGrid, x0, y0, R2x, R2y, xradius, yradius)
    xi = findall(x -> abs(x - x0) ≤ xradius, g.x)
    yi = findall(y -> abs(y - y0) ≤ yradius, g.y)
    zx = NMRTools.NMRBase._lineshape(2π * hz(x0, g.xaxis), R2x, g.ωx[xi],
                                     g.xaxis[:window], RealLineshape())
    zy = NMRTools.NMRBase._lineshape(2π * hz(y0, g.yaxis), R2y, g.ωy[yi],
                                     g.yaxis[:window], RealLineshape())
    return xi, yi, (π^2 * R2x * R2y) .* (zx .* zy')
end

"""
    samegrid(expt, i, j) -> Bool

Whether planes `i` and `j` share chemical shifts, frequencies and apodisation, so that a
peak has the same shape in both.
"""
function samegrid(expt::Experiment, i, j)
    s = expt.specdata
    s.nmrdata[i] === s.nmrdata[j] && return true
    s.x[i] == s.x[j] && s.y[i] == s.y[j] || return false
    for dim in (F1Dim, F2Dim)
        a = dims(s.nmrdata[i], dim)
        b = dims(s.nmrdata[j], dim)
        a[:window] == b[:window] || return false
        δ = data(s.nmrdata[i], dim)
        hz(first(δ), a) == hz(first(δ), b) && hz(last(δ), a) == hz(last(δ), b) ||
            return false
    end
    return true
end

"Position and linewidths of `peak` in plane `i`, as `(x, y, R2x, R2y)`."
function shapeparams(peak::Peak, i, quantity=:value)
    return map(POSITION_PARAMS) do sym
        par = peak.parameters[sym]
        return (quantity === :initial ? par.initialvalue[] : par.value[])[i]
    end
end

"""
    simulate!(z, peak, expt)

Add `peak`, at its fitted amplitudes, to every plane of `z`. Its shape is computed once
for each run of planes sharing a grid, rather than once per plane.
"""
function simulate!(z, peak::Peak, expt::Experiment)
    amp = peak.parameters[:amp].value[]
    cached = nothing
    for i in 1:nslices(expt)
        if isnothing(cached) || !hasfixedpositions(expt) || !samegrid(expt, i, i - 1)
            radii = (peak.xradius[], peak.yradius[])
            cached = peakshape(PlaneGrid(expt, i), shapeparams(peak, i)..., radii...)
        end
        xi, yi, block = cached
        z[i][xi, yi] .+= amp[i] .* block
    end
    return z
end

# ---- the cluster fit --------------------------------------------------------------

"""
    FitGroup

Planes of a cluster sharing one grid and one mask, and so one design matrix, whose column
`p` is peak `p` at unit amplitude over the masked points. `z` holds the observed points of
each plane, and `active` marks the planes that constrain the lineshapes; a skipped plane
still has its amplitudes measured, at the shapes the others determined.
"""
struct FitGroup
    grid::PlaneGrid
    mask::BitVector
    planes::Vector{Int}
    z::Matrix{Float64}
    active::BitVector
end

"""
    fitgroups(expt, masks, planes; skip=Set{Int}()) -> Vector{FitGroup}

Gather `planes` into groups sharing a grid and a mask, given each plane's cluster mask.
Every plane is active if `skip` would otherwise leave none.
"""
function fitgroups(expt::Experiment, masks, planes; skip=Set{Int}())
    all(in(skip), planes) && (skip = Set{Int}())
    reps = Int[]
    members = Vector{Int}[]
    for i in planes
        j = findfirst(r -> (masks[i] === masks[r] || masks[i] == masks[r]) &&
                           samegrid(expt, i, r), reps)
        if isnothing(j)
            push!(reps, i)
            push!(members, [i])
        else
            push!(members[j], i)
        end
    end
    return map(reps, members) do r, ps
        xb = vec(any(masks[r]; dims=2))
        yb = vec(any(masks[r]; dims=1))
        mask = BitVector(vec(masks[r][xb, yb]))
        z = stack(vec(expt.specdata.z[i][xb, yb])[mask] for i in ps)
        return FitGroup(PlaneGrid(expt, r, xb, yb), mask, ps, z,
                        BitVector([i ∉ skip for i in ps]))
    end
end

"""
    designmatrix(group, φ, radii) -> Matrix

The design matrix of `group` for peak shapes `φ` (x, y, R2x and R2y of each peak in turn).
"""
function designmatrix(g::FitGroup, φ, radii)
    A = zeros(count(g.mask), length(radii))
    block = zeros(length(g.grid.x), length(g.grid.y))
    for p in eachindex(radii)
        xi, yi, b = peakshape(g.grid, φ[4p - 3], φ[4p - 2], φ[4p - 1], φ[4p], radii[p]...)
        fill!(block, 0.0)
        block[xi, yi] .= b
        A[:, p] .= view(vec(block), g.mask)
    end
    return A
end

"Least-squares amplitudes, one column per plane, tolerating coincident peaks."
amplitudes(A, Z) = qr(A, ColumnNorm()) \ Z

"""
    ShapeFit

Peak shapes in fit coordinates. Positions are fitted as offsets from their starting point in
units of the fitting radius, bounded by ±1, and linewidths as they are. A chemical shift of
order a hundred ppm would otherwise set the scale of the relative convergence tolerance.
"""
struct ShapeFit
    origin::Vector{Float64}
    scale::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

physical(f::ShapeFit, θ) = f.origin .+ f.scale .* θ

"""
    maxR2(axis, δ, radius) -> Float64

The broadest linewidth the fitting window can determine: a Lorentzian whose full width at
half height, R2/π, spans twice the window (four radii). Anything broader is barely curved
across the window, so widening the radius is what allows a broader line.
"""
maxR2(axis, δ, radius) = 4π * abs(hz(δ + radius, axis) - hz(δ, axis))

function ShapeFit(peaks, i, grid::PlaneGrid)
    origin = Float64[]
    scale = Float64[]
    lower = Float64[]
    upper = Float64[]
    for peak in peaks
        x0, y0, _, _ = shapeparams(peak, i, :initial)
        append!(origin, (x0, y0, 0.0, 0.0))
        append!(scale, (peak.xradius[], peak.yradius[], 1.0, 1.0))
        append!(lower,
                (-1.0, -1.0, peak.parameters[:R2x].minvalue[],
                 peak.parameters[:R2y].minvalue[]))
        append!(upper,
                (1.0, 1.0,
                 max(peak.parameters[:R2x].maxvalue[],
                     maxR2(grid.xaxis, x0, peak.xradius[])),
                 max(peak.parameters[:R2y].maxvalue[],
                     maxR2(grid.yaxis, y0, peak.yradius[]))))
    end
    return ShapeFit(origin, scale, lower, upper)
end

"""
    fitshapes(peaks, groups, i, check; warm=false) -> NamedTuple

Fit the shapes of `peaks` in plane `i` (the shared shapes of a fixed-peak series, or one
plane of a moving-peak series) to `groups`, returning the shapes `φ` and their errors, each
plane's amplitudes and errors, and the fit's status. A cold start begins from the placed
positions and the initial linewidths, a warm one from the last fitted values. `check` is
called at every step and throws to stop the fit; nothing is written to the peaks here.
"""
function fitshapes(peaks, groups, i, check; warm=false)
    f = ShapeFit(peaks, i, first(groups).grid)
    radii = [(p.xradius[], p.yradius[]) for p in peaks]
    φ0 = reduce(vcat, [collect(shapeparams(p, i, warm ? :value : :initial)) for p in peaks])
    θ0 = clamp.((φ0 .- f.origin) ./ f.scale, f.lower, f.upper)
    fitted = [(g, g.z[:, g.active]) for g in groups if any(g.active)]

    function resid(θ)
        check()
        φ = physical(f, θ)
        r = Float64[]
        for (g, Z) in fitted
            A = designmatrix(g, φ, radii)
            append!(r, vec(Z .- A * amplitudes(A, Z)))
        end
        return r
    end

    sol = LsqFit.lmfit(resid, θ0, Float64[]; lower=f.lower, upper=f.upper,
                       autodiff=:finite, maxIter=200, x_tol=1e-8, g_tol=1e-10)
    θ = coef(sol)
    check()
    φ = physical(f, θ)
    σθ, amps, σamps = uncertainties(groups, f, θ, radii)

    atbound = any(j -> θ[j] ≤ f.lower[j] + 1e-6 * abs(f.lower[j]) + 1e-9 ||
                       θ[j] ≥ f.upper[j] - 1e-6 * abs(f.upper[j]) - 1e-9,
                  eachindex(θ))
    status = !sol.converged ? :maxiter : atbound ? :bound : :converged
    return (; φ, σφ=f.scale .* σθ, amps, σamps, status)
end

"""
    uncertainties(groups, f, θ, radii) -> (σθ, amps, σamps)

Standard errors of the shapes and of every plane's amplitudes, from the covariance of the
full problem, σ²(JᵀJ)⁻¹ with σ² estimated from the residuals. The amplitude blocks of JᵀJ
are small and independent, so the shape covariance is the inverse of their Schur
complement, S = Σᵢ Dᵢᵀ(I − A(AᵀA)⁻¹Aᵀ)Dᵢ over the active planes, where Dᵢ is the derivative
of A·aᵢ with respect to the shapes. Each amplitude's covariance is then
σ²[(AᵀA)⁻¹ + (AᵀA)⁻¹AᵀDᵢS⁻¹DᵢᵀA(AᵀA)⁻¹], its own noise plus that carried from the shapes,
which also holds for a skipped plane. `amps` and `σamps` map each plane to one value per
peak. Errors are `NaN` where the data do not determine them.
"""
function uncertainties(groups, f::ShapeFit, θ, radii)
    nθ = length(θ)
    φ = physical(f, θ)
    S = zeros(nθ, nθ)
    rss = 0.0
    npoints = 0
    nlinear = 0
    cached = map(groups) do g
        A = designmatrix(g, φ, radii)
        a = amplitudes(A, g.z)
        # dA[j] is the derivative of the one column of A that shape parameter j changes
        dA = map(1:nθ) do j
            h = 1e-6 * max(abs(θ[j]), 1.0)
            θp = copy(θ)
            θm = copy(θ)
            θp[j] += h
            θm[j] -= h
            p = cld(j, 4)
            return (designmatrix(g, physical(f, θp), radii)[:, p] .-
                    designmatrix(g, physical(f, θm), radii)[:, p]) ./ (2h)
        end
        Ginv = covinverse(A' * A)
        Ds = map(eachindex(g.planes)) do k
            return stack(dA[j] .* a[cld(j, 4), k] for j in 1:nθ)
        end
        for k in findall(g.active)
            D = Ds[k]
            S .+= D' * (D .- A * (Ginv * (A' * D)))
            rss += sum(abs2, g.z[:, k] .- A * a[:, k])
            npoints += size(g.z, 1)
            nlinear += size(A, 2)
        end
        return (A, a, Ginv, Ds)
    end
    mse = rss / max(npoints - nθ - nlinear, 1)
    Sinv = covinverse(S)
    σθ = sqrt.(max.(mse .* diag(Sinv), 0.0))

    amps = Dict{Int,Vector{Float64}}()
    σamps = Dict{Int,Vector{Float64}}()
    for (g, (A, a, Ginv, Ds)) in zip(groups, cached)
        for (k, i) in enumerate(g.planes)
            B = Ginv * (A' * Ds[k])
            C = Ginv .+ B * Sinv * B'
            amps[i] = a[:, k]
            σamps[i] = sqrt.(max.(mse .* diag(C), 0.0))
        end
    end
    return σθ, amps, σamps
end

"Inverse of a symmetric positive-definite matrix, or `NaN`s where it is singular."
function covinverse(M)
    F = cholesky(Symmetric(M); check=false)
    return issuccess(F) ? inv(F) : fill(NaN, size(M))
end

"Write a shape fit into `peaks`: their shapes in plane `i`, and their amplitudes."
function commit!(peaks, i, result)
    for (p, peak) in enumerate(peaks)
        for (j, sym) in enumerate(POSITION_PARAMS)
            peak.parameters[sym].value[][i] = result.φ[4(p - 1) + j]
            peak.parameters[sym].uncertainty[][i] = result.σφ[4(p - 1) + j]
        end
        for (plane, a) in result.amps
            peak.parameters[:amp].value[][plane] = a[p]
            peak.parameters[:amp].uncertainty[][plane] = result.σamps[plane][p]
        end
    end
    return peaks
end

# Statuses ranked so that a moving-peak cluster reports the worst of its planes.
const STATUS_RANK = Dict(:unfitted => 0, :converged => 1, :bound => 2, :maxiter => 3,
                         :timeout => 4)
worststatus(a, b) = STATUS_RANK[a] ≥ STATUS_RANK[b] ? a : b

"""
    fitcluster!(peaks, expt, check; warm=false) -> Symbol

Fit one cluster of overlapping peaks and write the results into them, returning how the fit
ended. A fixed-peak series is fitted jointly across its planes, sharing one set of shapes; a
moving-peak series is fitted plane by plane, since its shapes change from one to the next.
"""
function fitcluster!(peaks, expt::Experiment, check; warm=false)
    masks = mask(peaks, expt)
    groups = fitgroups(expt, masks, 1:nslices(expt); skip=skipset(expt))
    result = fitshapes(peaks, groups, 1, check; warm)
    commit!(peaks, 1, result)
    return result.status
end

function fitcluster!(peaks, expt::MovingPeakExperiment, check; warm=false)
    masks = mask(peaks, expt)
    results = map(1:nslices(expt)) do i
        return fitshapes(peaks, fitgroups(expt, masks, [i]), i, check; warm)
    end
    # nothing is written until every plane has finished, so a stopped fit changes nothing
    for (i, result) in enumerate(results)
        commit!(peaks, i, result)
    end
    return reduce(worststatus, (r.status for r in results); init=:converged)
end

"""
    fitcheck(expt, generation, budget) -> Function

A check for a fit started at `generation`: it throws `FitCancelled` once the fit has been
superseded and `FitTimeout` once `budget` seconds have passed. Single-threaded, it also
yields so that the window can register a cancellation.
"""
function fitcheck(expt, generation, budget)
    t0 = time()
    return function ()
        expt.state[][:fit_generation][] == generation || throw(FitCancelled())
        time() - t0 > budget && throw(FitTimeout())
        Threads.nthreads() == 1 && yield()
        return nothing
    end
end

"""
    fitclusters!(expt, clusters, generation, budget; warm=false) -> Vector{Symbol}

Fit each cluster (a vector of peaks), in parallel when Julia has more than one thread, and
return how each ended, `:cancelled` included. Progress is published to
`state[:fitprogress]`.

Call this from the main thread. The worker threads only compute and write into the peaks'
parameter arrays; every observable, and so every plot, is updated from the calling task,
since GLMakie cannot update a plot from another thread.
"""
function fitclusters!(expt, clusters, generation, budget; warm=false)
    progress = expt.state[][:fitprogress]
    n = length(clusters)
    done = Threads.Atomic{Int}(0)
    function run(peaks)
        status = try
            fitcluster!(peaks, expt, fitcheck(expt, generation, budget); warm)
        catch e
            e isa FitTimeout ? :timeout : e isa FitCancelled ? :cancelled : rethrow()
        end
        Threads.atomic_add!(done, 1)
        return status
    end
    if Threads.nthreads() == 1
        return map(clusters) do peaks
            status = run(peaks)
            progress[] = (done[], n)
            return status
        end
    end
    tasks = [Threads.@spawn run(peaks) for peaks in clusters]
    while !all(istaskdone, tasks)
        progress[] = (done[], n)
        sleep(0.05)
    end
    progress[] = (done[], n)
    return fetch.(tasks)
end
