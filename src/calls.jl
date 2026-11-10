# The call that produced an analysis, recorded by each entry point and saved with its output
# as `analysis.toml`, which `reopen` reads to repeat it. Shared by every analysis module.

"""
    AnalysisCall

The call that produced an analysis: the entry point's name, its arguments and keywords as
resolved, and the working directory they were given in, against which relative paths
resolve. Written to `analysis.toml` when the analysis is saved; see [`reopen`](@ref).
"""
struct AnalysisCall
    func::String
    args::Vector{Any}
    kwargs::Vector{Pair{Symbol,Any}}
    dir::String
end

"""
    analysiscall(func, args...; kwargs...) -> AnalysisCall

Record a call. Keyword arguments whose values have no literal (a `SeriesModel` object, say)
are left out, since the point of the record is that it can be read back.
"""
function analysiscall(func::AbstractString, args...; kwargs...)
    kws = Pair{Symbol,Any}[k => v for (k, v) in pairs(kwargs) if !isnothing(callvalue(v))]
    return AnalysisCall(String(func), collect(Any, args), kws, pwd())
end

"""
`false` where a positional argument has no literal - a spectrum passed as an `NMRData`
rather than as a path - so that nothing misleading is written.
"""
reproducible(call::AnalysisCall) = all(a -> !isnothing(callvalue(a)), call.args)

"""
    callvalue(x) -> String or nothing

`x` as a Julia literal, or `nothing` where it has none. Vectors are written out in full: a
gradient ramp or a delay list is exactly what has to be repeated.
"""
callvalue(x::Union{Real,AbstractString,Symbol}) = repr(x)
callvalue(::Nothing) = nothing
callvalue(x) = nothing
function callvalue(x::Union{AbstractVector,Tuple})
    vs = map(callvalue, collect(x))
    any(isnothing, vs) && return nothing
    inner = join(vs, ", ")
    return x isa Tuple ? (length(x) == 1 ? "($inner,)" : "($inner)") : "[$inner]"
end
function callvalue(x::NamedTuple)
    vs = map(callvalue, values(x))
    any(isnothing, vs) && return nothing
    return "(" * join(("$k=$v" for (k, v) in zip(keys(x), vs)), ", ") *
           (length(x) == 1 ? ",)" : ")")
end
function callvalue(x::Pair)
    a, b = callvalue(x.first), callvalue(x.second)
    return isnothing(a) || isnothing(b) ? nothing : "$a => $b"
end

"""
    relocate(x, from, to; data=false)

`x` with every string naming a file or folder (relative to `from`) rewritten relative to
the folder `to`: as `./name` inside it, `../name` outside it. With `data`, which marks a
positional argument, every string is a path and is rewritten, even one that names no single
file, such as the nmrPipe template `ft/test%03d.ft2`; an experiment number naming a folder
is rewritten too.
"""
relocate(x, from, to; data=false) = x
function relocate(x::AbstractString, from, to; data=false)
    p = joinpath(from, x)
    data || ispath(p) || return x
    r = relpath(p, to)
    return startswith(r, "..") ? r : "./" * r
end
function relocate(x::Integer, from, to; data=false)
    return data && ispath(joinpath(from, string(x))) ? relocate(string(x), from, to) : x
end
relocate(x::Bool, from, to; data=false) = x
relocate(x::AbstractVector, from, to; data=false) = map(y -> relocate(y, from, to; data), x)
relocate(x::Tuple, from, to; data=false) = map(y -> relocate(y, from, to; data), x)

"""
    callstring(call, folder=nothing; extra...) -> String or nothing

`call` as Julia code, one keyword to a line, with the keywords `extra` appended. Given a
`folder`, paths are written relative to it, as they are in `analysis.toml`. `nothing` where
the call has no literal form (see [`reproducible`](@ref)).
"""
function callstring(call::AnalysisCall, folder=nothing; extra...)
    reproducible(call) || return nothing
    move(x; data=false) = isnothing(folder) ? x : relocate(x, call.dir, folder; data)
    args = [callvalue(move(a; data=true)) for a in call.args]
    kwargs = [k => callvalue(move(v)) for (k, v) in [call.kwargs; collect(pairs(extra))]]
    filter!(kv -> !isnothing(kv.second), kwargs)
    io = IOBuffer()
    print(io, call.func, "(", join(args, ", "))
    if isempty(kwargs)
        print(io, ")")
    else
        isempty(args) || print(io, ";")
        println(io)
        pad = " "^(length(call.func) + 1)
        for (i, (name, value)) in enumerate(kwargs)
            print(io, pad, name, "=", value, i == length(kwargs) ? ")" : ",\n")
        end
    end
    return String(take!(io))
end

"""
The output folder a window starts with: `out` in the working directory, or, while
[`reopen`](@ref) repeats an analysis, the folder it reopened.
"""
const OUTPUT = Ref("out")
outputfolder() = OUTPUT[]

# An example's saved analysis is inside a downloaded artifact, which must not be changed.
function inartifact(path)
    return any(d -> startswith(abspath(path), joinpath(abspath(d), "artifacts")),
               DEPOT_PATH)
end

"""The file an output folder records its analysis in, which [`reopen`](@ref) reads."""
const RECORD = "analysis.toml"

"""
    writereproduce(io, folder, call; extra...)

Save the record of `call` in `folder` as `analysis.toml`, with the keywords `extra` (such
as the saved `peaklist`) added, and print the `Reproduce:` block of `summary.txt` to `io`.
Nothing is written for a call with no literal form, or for no call.
"""
writereproduce(io, folder, ::Nothing; extra...) = nothing
function writereproduce(io, folder, call::AnalysisCall; extra...)
    text = callstring(call, abspath(folder); extra...)
    isnothing(text) && return nothing
    open(joinpath(folder, RECORD), "w") do f
        println(f, "# The analysis that wrote this folder, with paths relative to it.")
        println(f, "# reopen(\"<this folder>\") in NMRAnalysis.jl repeats it.")
        println(f, "nmranalysis = ", repr(string(pkgversion(@__MODULE__))))
        println(f, "call = '''")
        println(f, text)
        return println(f, "'''")
    end
    println(io)
    println(io, "Reproduce:")
    println(io, "    reopen(", repr(relpath(abspath(folder))), ")")
    return nothing
end

"""
    reopen(path; kwargs...)

Reopen a saved analysis: repeat the call recorded in the output folder `path` (or the
`analysis.toml` file itself), with its parameters as they were resolved and the peak list
it saved, so the window opens as it was left and the fit repeats.

Keywords given here replace those recorded: `peaklist=nothing` reopens the same data with
the same parameters but no peaks. The window saves back into the folder reopened, unless it
is an example's saved analysis, which must not change: that is saved to `out`.

The record is read as data, never run as code: only literal values are accepted.

# Examples
```julia
reopen("out")                                    # an analysis saved earlier
reopen(exampleanalysis("relaxation2d_asyn"))     # an example's saved analysis
reopen("out"; peaklist=nothing)                  # the same, picking the peaks afresh
```
"""
function reopen(path::AbstractString; kwargs...)
    file = isdir(path) ? joinpath(path, RECORD) : path
    isfile(file) || throw(ArgumentError("no saved analysis at $path: no $RECORD"))
    folder = dirname(abspath(file))
    ex = Meta.parse(TOML.parsefile(file)["call"])
    (ex isa Expr && ex.head == :call && ex.args[1] isa Symbol) ||
        throw(ArgumentError("$file does not record a call"))
    f = ex.args[1]
    (f in names(@__MODULE__) && isdefined(@__MODULE__, f)) ||
        throw(ArgumentError("$file records $f, which is not an NMRAnalysis analysis"))
    args = Any[]
    kws = Pair{Symbol,Any}[]
    for a in ex.args[2:end]
        if a isa Expr && a.head == :parameters
            append!(kws, [kw.args[1] => literal(kw.args[2]) for kw in a.args])
        elseif a isa Expr && a.head == :kw
            push!(kws, a.args[1] => literal(a.args[2]))
        else
            push!(args, literal(a))
        end
    end
    inside(x) = resolve(x, folder)
    recorded = (; (k => inside(v) for (k, v) in kws)...)
    previous = OUTPUT[]
    inartifact(folder) || (OUTPUT[] = relpath(folder))
    try
        return getfield(@__MODULE__, f)(map(inside, args)...; merge(recorded, kwargs)...)
    finally
        OUTPUT[] = previous
    end
end

"""
A recorded path, written relative to the record (`./…`, `../…`), made absolute. Throws an
`ArgumentError` if it is not there, as when the folder was moved after it was saved.
"""
resolve(x, folder) = x
function resolve(x::AbstractString, folder)
    startswith(x, "./") || startswith(x, "../") || return x
    p = normpath(joinpath(folder, x))
    # an nmrPipe plane series is recorded as a template, such as ft/test%03d.ft2
    ispath(p) || occursin('%', basename(p)) ||
        throw(ArgumentError("$folder records $x, but there is nothing at $p. Was the " *
                            "folder moved, or the data renamed, after the analysis was saved?"))
    return p
end
resolve(x::AbstractVector, folder) = map(y -> resolve(y, folder), x)
resolve(x::Tuple, folder) = map(y -> resolve(y, folder), x)

"""
    literal(ex) -> value

The value of a literal expression: a number, string, symbol, `nothing`, or a vector, tuple,
named tuple or pair of literals. Anything else is refused, so a record is never run as code.
"""
literal(x::Union{Number,AbstractString}) = x
literal(x::QuoteNode) = x.value isa Symbol ? x.value : refuse(x)
function literal(x::Symbol)
    return x === :nothing ? nothing : x in (:Inf, :NaN) ? getfield(Base, x) : refuse(x)
end
function literal(ex::Expr)
    if ex.head == :vect
        return map(identity, Any[literal(a) for a in ex.args])
    elseif ex.head == :tuple && all(a -> a isa Expr && a.head == :(=), ex.args)
        return (; (a.args[1] => literal(a.args[2]) for a in ex.args)...)
    elseif ex.head == :tuple
        return Tuple(literal(a) for a in ex.args)
    elseif ex.head == :call && ex.args[1] == :(=>) && length(ex.args) == 3
        return literal(ex.args[2]) => literal(ex.args[3])
    elseif ex.head == :call && ex.args[1] == :- && length(ex.args) == 2
        v = literal(ex.args[2])
        v isa Number || refuse(ex)
        return -v
    end
    return refuse(ex)
end
literal(x) = refuse(x)

refuse(x) = throw(ArgumentError("not a literal value, so not read from a record: $x"))
