# Interactive parameter entry.
#
# Every 1D routine resolves a parameter the same way, in one fixed order of precedence:
#
#   1. an explicit argument, if the caller supplied one;
#   2. whatever the spectrum records - a pulse-sequence annotation first, then an
#      acquisition parameter;
#   3. failing both, a question asked here, before the analysis window opens.
#
# The entry points express that as `@something(explicit, annotation(...), acqusvalue(...),
# ask...(; prompt))`, which short-circuits, so a question is only ever printed for a value
# that genuinely could not be found. `prompt=false` (the default outside an interactive
# session) turns each question into either its stated default or an informative error, so a
# script or a test never blocks on stdin.
#
# Nothing here knows anything about NMR: the physics of *which* parameters an experiment
# needs, and what they mean, stays in that experiment's own `expt-*.jl`.

"""
    ask(label, default; unit="", note="", type=Float64, prompt=true) -> value
    ask(label; unit="", note="", type=Float64, prompt=true) -> value

Ask the user for a parameter value, on one line.

Given a `default`, the value is offered for confirmation and pressing enter accepts it;
without one (or with `default === nothing`) a value must be typed. `unit` and `note` are
shown alongside to say what is wanted and where the default came from. Re-asks until the
answer parses as `type`.

With `prompt=false` nothing is printed: the `default` is returned, or an `ArgumentError`
thrown if there is none.
"""
function ask(label::AbstractString, default; unit::AbstractString="",
             note::AbstractString="", type::Type=Float64, prompt::Bool=true)
    prompt || return default
    suffix = isempty(unit) ? "" : " $unit"
    while true
        print("$label = $default$suffix$note. Press enter to accept, or type a new value: ")
        response = strip(readline())
        isempty(response) && return default
        value = tryparse(type, response)
        isnothing(value) || return value
        println("  Sorry, \"$response\" is not a valid number. Please try again.")
    end
end

function ask(label::AbstractString; unit::AbstractString="", note::AbstractString="",
             type::Type=Float64, prompt::Bool=true)
    prompt ||
        throw(ArgumentError("$label could not be read from the spectrum - " *
                            "please supply it explicitly"))
    while true
        print("Enter $label$(isempty(unit) ? "" : " / $unit")$note: ")
        response = strip(readline())
        value = isempty(response) ? nothing : tryparse(type, response)
        isnothing(value) || return value
        println("  Sorry, that is not a valid number. Please try again.")
    end
end

# A missing default is the same question as no default at all - so an entry point can pass
# whatever it found (or didn't) straight through, without branching on it first.
function ask(label::AbstractString, ::Nothing; kwargs...)
    return ask(label; kwargs...)
end

"""
    askchoice(question, options; default=1, prompt=true) -> Int

Ask the user to pick one of `options`, returning its index. With `prompt=false`, or if the
menu is cancelled, `default` is returned.
"""
function askchoice(question::AbstractString, options; default::Integer=1,
                   prompt::Bool=true)
    prompt || return default
    println(question)
    choice = request(RadioMenu(collect(String, options)))
    return choice < 1 ? default : choice
end

"""
    askvector(label, n; unit="", prompt=true) -> Vector{Float64}

Ask for the `n` values of an arrayed acquisition variable - a list of relaxation delays,
say. They can be typed directly, separated by spaces or commas, or the name given of a
file holding one per line, which is how Bruker writes a `vdlist`.
"""
function askvector(label::AbstractString, n::Integer; unit::AbstractString="",
                   prompt::Bool=true)
    prompt ||
        throw(ArgumentError("$label could not be read from the spectrum - " *
                            "please supply them explicitly"))
    suffix = isempty(unit) ? "" : " (in $unit)"
    println("The experiment has $n spectra, but the $label are not recorded in the data.")
    while true
        print("Enter $n $label$suffix, separated by spaces, " *
              "or a file holding one per line: ")
        vals = parsevector(strip(readline()))
        if isnothing(vals)
            println("  Sorry, that could not be read as a list of numbers or as a file. " *
                    "Please try again.")
        elseif length(vals) != n
            println("  That is $(length(vals)) value(s), but the experiment has " *
                    "$n spectra. Please try again.")
        else
            return vals
        end
    end
end

"""
    parsevector(response) -> Vector{Float64} or nothing

Read a whitespace- or comma-separated list of numbers, or the contents of a file named by
`response`. Returns `nothing` if the response is neither.
"""
function parsevector(response::AbstractString)
    isempty(response) && return nothing
    isfile(response) && return parsevector(join(readlines(response), " "))
    fields = split(response, r"[\s,]+"; keepempty=false)
    vals = tryparse.(Float64, fields)
    return any(isnothing, vals) ? nothing : Vector{Float64}(vals)
end

"""
    askpath(label) -> String

Ask for the path of an experiment folder, re-asking until it exists. Used by the entry
points that can be called with no arguments at all.
"""
function askpath(label::AbstractString)
    while true
        print("Enter path to $label (i.e. a Bruker experiment folder): ")
        path = strip(readline())
        ispath(path) && return String(path)
        println("  No such file or directory: \"$path\". Please try again.")
    end
end
