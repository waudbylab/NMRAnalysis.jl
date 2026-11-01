# Shared conventions for the interactive analysis windows.

"""
    windowname(analysis) -> String

The title of an analysis window, naming the package and the analysis:
`windowname("1D peak tracking") == "NMRAnalysis.jl - 1D peak tracking"`.
"""
windowname(analysis) = "NMRAnalysis.jl - $analysis"

"A name as it reads mid-sentence: a lower-case first letter, unless it begins an acronym."
function midsentence(name)
    name = String(name)
    return length(name) > 1 && isuppercase(name[2]) ? name : lowercasefirst(name)
end
