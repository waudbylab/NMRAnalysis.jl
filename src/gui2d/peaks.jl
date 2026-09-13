function Peak(initialposition, label, xradius=0.03, yradius=0.3)
    xy = MaybeVector(initialposition)
    x = [pos[1] for pos in xy]
    y = [pos[2] for pos in xy]
    length(x)
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
                postpars,
                Observable(false)) # post-fitted
end

function position(peak::Peak)
    return lift((x, y) -> MaybeVector(Point2f.(x, y)), peak.parameters[:x].value,
                peak.parameters[:y].value)
end
function initialposition(peak::Peak)
    return lift((x, y) -> MaybeVector(Point2f.(x, y)), peak.parameters[:x].initialvalue,
                peak.parameters[:y].initialvalue)
end

# generic overlap function - can specialise for different experiments
isoverlapping(peak1, peak2, ::Experiment) = isoverlapping(peak1, peak2)
function isoverlapping(peak1, peak2)
    # Δδs will be a MaybeVector - could be a single element or a list of different shifts
    Δδs = MaybeVector(initialposition(peak1)[] .- initialposition(peak2)[])
    return any(Δδ -> begin
                   dX = Δδ[1] / (peak1.xradius[] + peak2.xradius[])
                   dY = Δδ[2] / (peak1.yradius[] + peak2.yradius[])
                   (dX^2 + dY^2) <= 1.0
               end, Δδs)
end

function pack(peaks::Vector{Peak}, quantity=:value)
    p = Vector{Float64}()
    for peak in peaks
        pack!(p, peak, quantity)
    end
    return p
end

"Create vector of parameters (options - :min, :max, :initial)"
function pack!(p, peak::Peak, quantity=:value)
    d = peak.parameters
    if quantity == :min
        d[:x].minvalue[] = d[:x].initialvalue[] .- peak.xradius[]
        d[:y].minvalue[] = d[:y].initialvalue[] .- peak.yradius[]
    elseif quantity == :max
        d[:x].maxvalue[] = d[:x].initialvalue[] .+ peak.xradius[]
        d[:y].maxvalue[] = d[:y].initialvalue[] .+ peak.yradius[]
    end
    for par in values(d)
        pack!(p, par, quantity)
    end

    return p
end

function unpack!(v, peaks::Vector{Peak}, quantity=:value)
    @debug "Unpacking peaks" maxlog = 10
    for peak in peaks
        unpack!(v, peak, quantity)
    end
end

"Unpack a parameter and pop from input vector. Quantity could also be :uncertainty"
function unpack!(v, peak::Peak, quantity=:value)
    @debug "Unpacking peak $(peak.label)" maxlog = 10
    d = peak.parameters
    for par in values(d)
        unpack!(v, par, quantity)
    end
end

function Base.show(io::IO, p::Peak)
    return print(io, "Peak($(p.label[]), position=$(position(p)[]))")
end

function Base.show(io::IO, mime::MIME"text/plain", p::Peak)
    println(io, "Peak: $(p.label[])")
    println(io, "  initial position: $(initialposition(p)[])")
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
