module MaybeVectorModule

export MaybeVector, SingleElementVector, StandardVector

"""
    MaybeVector{T} <: AbstractVector{T}

Abstract type representing a vector that may contain either a single element or multiple elements.
All indices return the same value for SingleElementVector.
"""
abstract type MaybeVector{T} <: AbstractVector{T} end

"""
    SingleElementVector{T} <: MaybeVector{T}

Vector type that returns its single element for any valid index.
"""
mutable struct SingleElementVector{T} <: MaybeVector{T}
    x::T
end

"""
    StandardVector{T} <: MaybeVector{T}

Standard vector implementation wrapping Vector{T}.
"""
struct StandardVector{T} <: MaybeVector{T}
    x::Vector{T}
end

# Constructors
MaybeVector(x::T) where {T} = SingleElementVector{T}(x)
MaybeVector(x::Vector{T}) where {T} = StandardVector{T}(x)

# SingleElementVector interface
Base.size(mv::SingleElementVector) = (1,)
Base.getindex(mv::SingleElementVector, i::Int) = mv.x
Base.setindex!(mv::SingleElementVector, v, i::Int) = (mv.x = v)
Base.length(mv::SingleElementVector) = 1
Base.IndexStyle(::Type{<:SingleElementVector}) = IndexLinear()
Base.similar(mv::SingleElementVector{T}) where {T} = SingleElementVector{T}(mv.x)
Base.copy(mv::SingleElementVector) = SingleElementVector(mv.x)

# StandardVector interface
Base.size(mv::StandardVector) = size(mv.x)
Base.getindex(mv::StandardVector, i::Int) = mv.x[i]
Base.setindex!(mv::StandardVector, v, i::Int) = (mv.x[i] = v)
Base.length(mv::StandardVector) = length(mv.x)
Base.IndexStyle(::Type{<:StandardVector}) = IndexLinear()
Base.similar(mv::StandardVector{T}) where {T} = StandardVector{T}(similar(mv.x))
Base.copy(mv::StandardVector) = StandardVector(copy(mv.x))

# Iteration protocol
Base.iterate(mv::SingleElementVector) = (mv.x, nothing)
Base.iterate(mv::SingleElementVector, state) = nothing
Base.iterate(mv::StandardVector) = iterate(mv.x)
Base.iterate(mv::StandardVector, state) = iterate(mv.x, state)

# Broadcasting
Base.broadcastable(mv::SingleElementVector) = Ref(mv.x)
Base.broadcastable(mv::StandardVector) = mv.x

# Conversions
function Base.convert(::Type{StandardVector{T}}, x::SingleElementVector{T}) where {T}
    return StandardVector([x.x])
end
function Base.convert(::Type{SingleElementVector{T}}, x::StandardVector{T}) where {T}
    return length(x.x) == 1 ? SingleElementVector(first(x.x)) :
           throw(ArgumentError("Cannot convert StandardVector with length != 1 to SingleElementVector"))
end
Base.convert(::Type{Vector{T}}, x::SingleElementVector{T}) where {T} = [x.x]
Base.convert(::Type{Vector{T}}, x::StandardVector{T}) where {T} = x.x
function Base.convert(::Type{MaybeVector{T}}, x::Vector{T}) where {T}
    return length(x) == 1 ? SingleElementVector(first(x)) : StandardVector(x)
end

end
