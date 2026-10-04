# The call that produced an analysis, recorded by each entry point so that `summary.txt` can
# print a line repeating it. Shared by Analysis1D and GUI2D.

"""
    AnalysisCall

The call that produced an analysis, recorded by the entry point so `summary.txt` can print
a line that repeats it. `reproducible` is false where one of the positional arguments has
no Julia literal - a spectrum passed as an `NMRData` rather than as a path - in which case
no such line is printed at all rather than a misleading one.
"""
struct AnalysisCall
    func::String
    args::Vector{String}
    kwargs::Vector{Pair{Symbol,String}}
    reproducible::Bool
end

"""
    analysiscall(func, args...; kwargs...) -> AnalysisCall

Record a call. Keyword arguments whose values have no literal (a `SeriesModel` object, say)
are left out, since the point of the record is that it can be pasted back into Julia.
"""
function analysiscall(func::AbstractString, args...; kwargs...)
    rendered = Pair{Symbol,String}[]
    for (name, value) in pairs(kwargs)
        text = callvalue(value)
        isnothing(text) || push!(rendered, name => text)
    end
    args_ = [callvalue(a) for a in args]
    return AnalysisCall(String(func), String[something(a, "spec") for a in args_],
                        rendered, !any(isnothing, args_))
end

"""
    callvalue(x) -> String or nothing

`x` as a Julia literal, or `nothing` where it has none. Vectors of numbers are written out
in full: a gradient ramp or a delay list is exactly what has to be repeated, so abbreviating
it would defeat the purpose.
"""
callvalue(x::Union{Real,Bool,Symbol,AbstractString}) = repr(x)
callvalue(x::AbstractVector{<:Real}) = repr(collect(x))
# the experiment folders of a multi-file analysis: the one argument worth repeating in full
callvalue(x::AbstractVector{<:AbstractString}) = repr(collect(x))
callvalue(x::NamedTuple) = repr(x)
callvalue(::Nothing) = nothing
callvalue(x) = nothing

# A tuple (2D dimension weights, component names) and a list of name => value pairs (a custom
# model's starting values) are written as they are
callvalue(x::Tuple) = repr(x)
callvalue(x::AbstractVector{<:Pair}) = repr(collect(x))

"""
    callstring(call, extra=[]) -> String or nothing

`call` as Julia code to paste back, with the keyword arguments in `extra` (names and
already-rendered values) appended, one keyword to a line. `nothing` where the call has no
literal form (see [`AnalysisCall`](@ref)).
"""
function callstring(call::AnalysisCall,
                    extra::Vector{Pair{Symbol,String}}=Pair{Symbol,String}[])
    call.reproducible || return nothing
    kwargs = [call.kwargs; extra]
    io = IOBuffer()
    print(io, call.func, "(", join(call.args, ", "))
    if isempty(kwargs)
        print(io, ")")
    else
        isempty(call.args) || print(io, ";")
        println(io)
        pad = " "^(length(call.func) + 1)
        for (i, (name, value)) in enumerate(kwargs)
            print(io, pad, name, "=", value, i == length(kwargs) ? ")" : ",\n")
        end
    end
    return String(take!(io))
end
