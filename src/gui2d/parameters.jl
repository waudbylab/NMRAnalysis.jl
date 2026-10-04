function Parameter(label, initialvalue; minvalue=(-Inf), maxvalue=Inf, uncertainty=Inf)
    if !isa(initialvalue, MaybeVector)
        initialvalue = MaybeVector(initialvalue)
    end
    value = deepcopy(initialvalue)

    # ensure uncertainty matches size of initialvalue
    if length(uncertainty) != length(value)
        uncertainty = fill(uncertainty, length(value))
    end
    uncertainty = MaybeVector(uncertainty)

    # min/max are stored untyped so they can hold either a scalar bound or one per plane.
    # Positions are bounded by the fitting radius instead (see ShapeFit).
    return Parameter(label, Observable(value), Observable(uncertainty),
                     Observable(initialvalue), Observable{Any}(minvalue),
                     Observable{Any}(maxvalue))
end

function Base.show(io::IO, p::Parameter)
    return print(io,
                 "Parameter($(p.label), value=$(p.value[]), uncertainty=$(p.uncertainty[]))")
end

function Base.show(io::IO, mime::MIME"text/plain", p::Parameter)
    println(io, "Parameter: $(p.label)")
    println(io, "  value: $(p.value[])")
    println(io, "  uncertainty: $(p.uncertainty[])")
    println(io, "  initial value: $(p.initialvalue[])")
    return println(io, "  bounds: [$(p.minvalue[]), $(p.maxvalue[])]")
end
