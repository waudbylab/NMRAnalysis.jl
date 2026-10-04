function Peak(initialposition, label, xradius=0.03, yradius=0.3)
    xy = MaybeVector(initialposition)
    # Float64: stored as the Float32 of a Point2f, a fitted position written back here would
    # be rounded to about 1e-5 ppm
    x = Float64[pos[1] for pos in xy]
    y = Float64[pos[2] for pos in xy]
    if length(x) == 1
        x = MaybeVector(x[1])
        y = MaybeVector(y[1])
    end
    pars = OrderedDict{Symbol,Parameter}()
    postpars = OrderedDict{Symbol,Parameter}()
    pars[:x] = Parameter("δx", x)
    pars[:y] = Parameter("δy", y)
    return Peak(pars,
                Observable(label),
                Observable(true), # touched
                Observable(xradius), # xradius
                Observable(yradius), # yradius
                Observable(false), # customradius
                postpars,
                Observable(false), # post-fitted
                Observable(:unfitted)) # fitstatus
end

"""
    fitstatus(peak) -> Symbol

How the last lineshape fit of `peak` ended: `:unfitted`, `:converged`, `:maxiter` (stopped
at the iteration limit), `:bound` (a position or linewidth finished on its limit) or
`:timeout` (stopped at the time limit, leaving the previous values in place).
"""
fitstatus(peak::Peak) = peak.fitstatus[]

"Whether the last fit of `peak` ended anywhere other than a converged optimum."
isunfinished(peak::Peak) = fitstatus(peak) in (:maxiter, :bound, :timeout)

"Whether the last fit of `peak` was stopped before it finished, so continuing it helps."
iscontinuable(peak::Peak) = fitstatus(peak) in (:maxiter, :timeout)

"""
    setradius!(peak, xradius, yradius)

Give `peak` radii of its own, which changing the experiment's default radii then leaves
alone.
"""
function setradius!(peak::Peak, xradius, yradius)
    peak.xradius.val = xradius
    peak.yradius.val = yradius
    peak.customradius.val = true
    peak.touched.val = true
    return peak
end

"Fitted position of `peak`, one `Point2f` for every plane or one per plane."
function position(peak::Peak)
    return MaybeVector(Point2f.(peak.parameters[:x].value[], peak.parameters[:y].value[]))
end

"Where `peak` was placed, one `Point2f` for every plane or one per plane."
function initialposition(peak::Peak)
    return MaybeVector(Point2f.(peak.parameters[:x].initialvalue[],
                                peak.parameters[:y].initialvalue[]))
end

# generic overlap function - can specialise for different experiments
isoverlapping(peak1, peak2, ::Experiment) = isoverlapping(peak1, peak2)
function isoverlapping(peak1, peak2)
    # Δδs will be a MaybeVector - could be a single element or a list of different shifts
    Δδs = MaybeVector(initialposition(peak1) .- initialposition(peak2))
    return any(Δδ -> begin
                   dX = Δδ[1] / (peak1.xradius[] + peak2.xradius[])
                   dY = Δδ[2] / (peak1.yradius[] + peak2.yradius[])
                   (dX^2 + dY^2) <= 1.0
               end, Δδs)
end

function Base.show(io::IO, p::Peak)
    return print(io, "Peak($(p.label[]), position=$(position(p)))")
end

function Base.show(io::IO, mime::MIME"text/plain", p::Peak)
    println(io, "Peak: $(p.label[])")
    println(io, "  initial position: $(initialposition(p))")
    println(io, "  radius: [$(p.xradius[]), $(p.yradius[])]")
    println(io, "  parameters:")
    for (k, v) in p.parameters
        println(io, "    $k: $(v.value[])")
    end
    if !isempty(p.postparameters)
        println(io, "  post-fit parameters:")
        for (k, v) in p.postparameters
            println(io, "    $k: $(v.value[])")
        end
    end
end
