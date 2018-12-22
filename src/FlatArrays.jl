module FlatArrays

# Data structures are not exported because their constructors should not be
# directly called.
export
    flatten,
    isflat

using Base:
    @propagate_inbounds,
    OneTo

import Base:
    IndexStyle,
    IndexLinear,
    axes,
    eltype,
    length,
    ndims,
    size,
    firstindex,
    lastindex,
    getindex,
    setindex!,
    iterate,
    has_offset_axes,
    convert,
    unsafe_convert

macro _inline_meta(); Expr(:meta, :inline); end

"""

All flat array types are sub-types of `AbstractFlatArray{T,N}`.

"""
abstract type AbstractFlatArray{T,N} <: DenseArray{T,N} end

"""

`FlatVector` type stores a *flat* vector.  Never directly call the
constructor but rather `flatten(A,n)` to build an instance of `FlatVector`.

"""
struct FlatVector{T,A<:AbstractArray{T}} <: AbstractFlatArray{T,1}
    parent::A
    nelem::Int
end

"""

`FlatMatrix` type stores a *flat* matrix.  Never directly call the
constructor but rather `flatten(A,m,n)` to build an instance of
`FlatMatrix`.

"""
struct FlatMatrix{T,A<:AbstractArray{T}} <: AbstractFlatArray{T,2}
    parent::A
    nelem::Int
    dims::NTuple{2,Int}
end

# Extend base methods so that flat arrays behave like arrays.

eltype(::AbstractFlatArray{T,N}) where {T,N} = T
ndims(::AbstractFlatArray{T,N}) where {T,N} = N
length(A::AbstractFlatArray) = A.nelem
has_offset_axes(::AbstractFlatArray) = false
firstindex(::AbstractFlatArray) = 1
firstindex(::AbstractFlatArray{T,N}, d::Integer) where {T,N} =
    (@_inline_meta; (d % UInt) - 1 < N ? 1 : error("dimension out of range"))
lastindex(A::AbstractFlatArray) = length(A)
lastindex(A::AbstractFlatArray, d) = size(A, d)
IndexStyle(::Type{<:AbstractFlatArray}) = IndexLinear()

unsafe_convert(::Type{Ptr{Cvoid}}, A::AbstractFlatArray) =
    unsafe_convert(Ptr{Cvoid}, A.parent)
unsafe_convert(::Type{Ptr{T}}, A::AbstractFlatArray{T}) where {T} =
    unsafe_convert(Ptr{T}, A.parent)

# The following is not needed:
# Base.IteratorSize(::Type{<:AbstractFlatArray{<:Any,N}}) where {N} =
#     Base.HasShape{N}()
# Base.IteratorEltype(::Type{<:AbstractFlatArray{<:Any,N}}) where {N} =
#     Base.HasEltype()

iterate(A::AbstractFlatArray, i=1) =
    (@_inline_meta; (i % UInt) - 1 < length(A) ?
     (@inbounds A[i], i + 1) : nothing)

size(A::FlatVector) = (A.nelem,)
size(A::FlatMatrix) = A.dims
size(A::FlatVector, d::Integer) = (d == 1 ? A.nelem :
                                   d > 1 ? 1 :
                                   error("dimension out of range"))
size(A::FlatMatrix, d::Integer) = (d > 2 ? 1 :
                                   d > 0 ? A.dims[d] :
                                   error("dimension out of range"))

axes(A::FlatVector) = (OneTo(A.nelem),)
axes(A::FlatMatrix) = map(OneTo, A.dims)
axes(A::FlatVector, d::Integer) = (d == 1 ? OneTo(A.nelem) :
                                   d > 1 ? OneTo(1) :
                                   error("dimension out of range"))
axes(A::FlatMatrix, d::Integer) = (d > 2 ? OneTo(1) :
                                   d > 0 ? OneTo(A.dims[d]) :
                                   error("dimension out of range"))

strides(A::FlatVector) = (1,)
strides(A::FlatMatrix) = (1, size(A, 1))
stride(A::FlatVector, d::Integer) = (d == 1 ? 1 :
                                     d >= 2 ? length(A) :
                                     error("dimension out of range"))
stride(A::FlatMatrix, d::Integer) = (d == 1 ? 1 :
                                     d == 2 ? size(A, 1) :
                                     d >= 3 ? length(A) :
                                     error("dimension out of range"))

@inline @propagate_inbounds getindex(A::AbstractFlatArray, i) =
    getindex(A, _index(A, i))
@inline @propagate_inbounds getindex(A::AbstractFlatArray, i, j) =
    getindex(A, _index(A, i, j))
@inline @propagate_inbounds function getindex(A::AbstractFlatArray, i::Int)
    @boundscheck checkbounds(A.parent, i)
    @inbounds val = A.parent[i]
    return val
end

@inline @propagate_inbounds setindex!(A::AbstractFlatArray, val, i) =
    setindex!(A, val, _index(A, i))
@inline @propagate_inbounds setindex!(A::AbstractFlatArray, val, i, j) =
    setindex!(A, val, _index(A, i, j))
@inline @propagate_inbounds function setindex!(A::AbstractFlatArray, val, i::Int)
    @boundscheck checkbounds(A.parent, i)
    @inbounds A.parent[i] = val
    return val
end

@inline _index(A::AbstractFlatArray, i::Integer) = Int(i)
@inline _index(A::FlatMatrix, i::Int, j::Int) = i + A.dims[1]*(j - 1)
@inline _index(A::FlatMatrix, i::Integer, j::Integer) = _index(A, Int(i), Int(j))
@inline _index(A::FlatVector, I::CartesianIndex{1}) = I[1]
@inline _index(A::FlatMatrix, I::CartesianIndex{2}) = _index(A, I[1], I[2])

"""

```julia
flatten(A, n) -> V
```

yields a *flat* vector `V` of length `n` whose elements are given by array `A`.

```julia
flatten(A, m, n) -> M
```

yields a *flat* matrix `M` with `m` rows and `n` columns whose elements are
given by array `A`.

A *flat* array has 1-based linear indexing, column-major storage order and all
its elements are contiguous.  This may be very convenient for calling efficient
code which expects such kind of arrays.

If the argument `A` of `flatten` is already a *flat* array, the returned array
will share the contents of `A`; otherwise an independent copy of the elements
is made.

If the contents of `A` can be shared, calling `flatten` is faster (by a factor
~ 7 for a vector, ~ 5 for a matrix) and uses less memory than calling
[`reshape`](@ref).  This may be significant for very small arrays.  Also most
base methods like `firstindex`, `lastindex`, etc. applied to an
`AbstractFlatArray` can be optimized out at compilation time.  Finally, the
array returned by `flatten` is guaranteed to have 1-based linear indexing and
column-major storage order.  However, contrary to `resize`, the array returned
by `flatten` may not share its contents with `A`.  In fact, if `isflat(A)` is
true, then `flatten(A,...)` will return an array sharing its contents with `A`;
otherwise, the contents is not shared.

!!! warning
    Many fast operations on flat arrays rely on the fact that the parent array
    backing the storage of the flat array is never resized (*e.g.*, with the
    [`resize!`](@ref) method).

See also: [`isflat`](@ref), [`resize`](@ref).

""" flatten

flatten(arr::AbstractArray, dims::NTuple{1,Integer}) =
    flatten(arr, dims[1])
flatten(arr::AbstractArray, nelem::Integer) =
    flatten(arr, Int(nelem))
@inline flatten(arr::A, nelem::Int) where {T,A<:AbstractArray{T}} =
    FlatVector(_flatten(arr, nelem), nelem)

flatten(arr::AbstractArray, dims::NTuple{2,Integer}) =
    flatten(arr, dims[1], dims[2])
flatten(arr::AbstractArray, nrows::Integer, ncols::Integer) =
    flatten(arr, Int(nrows), Int(ncols))
@inline function flatten(arr::A,
                         nrows::Int, ncols::Int) where {T,A<:AbstractArray{T}}
    nrows > 0 && ncols > 0 || error("invalid dimensions")
    nelem = nrows*ncols
    return FlatMatrix(_flatten(arr, nelem), nelem, (nrows, ncols))
end

@inline function _flatten(A::Array, nelem::Int)
    (nelem % UInt) - 1 < length(A) || error("incompatible length")
    return A
end

@inline function _flatten(A::AbstractFlatArray, nelem::Int)
    (nelem % UInt) - 1 < length(A) || error("incompatible length")
    return A.parent
end

@inline function _flatten(A::AbstractArray, nelem::Int)
    (nelem % UInt) - 1 < length(A) || error("incompatible length")
    return (isflat(A) ? A : _flatcopy(A, nelem))
end

function _flatcopy(A::AbstractArray{T}, nelem::Int) :: Vector{T} where {T}
    F = Array{T,1}(undef, nelem)
    F[1], state = iterate(A)
    @inbounds @simd for j in 2:nelem
        F[j], state = iterate(A, state)
    end
    return F
end

"""

`isflat(A)` yields whether array `A` is a *flat* array, that is an array
with 1-based indices and whose elements are all continuous.

See also: [`flatten`](@ref).

"""
isflat(::AbstractFlatArray) = true
isflat(::Array) = true
isflat(A::T) where {T<:AbstractArray} = _isflat(IndexStyle(T), A)

# The test below could also be `has_offset_axes(A) == false`
_isflat(::IndexLinear, A::AbstractArray) = (firstindex(A) == 1)

function _isflat(::IndexCartesian, A::AbstractArray{T,N}) where {T,N}
    #
    # We could have checked the strides with the following code:
    #
    #     n = 1
    #     @inbounds for d in 1:N
    #         stride(A, d) == n || return false
    #         n *= size(A, d)
    #     end
    #     return true
    #
    # However, (i) not all array types implement the stride method, (ii)
    # non-linear indexing style indicates that the array is certainly
    # non-flat (otherwise indexation is not optimized for that kind of array).
    # Therefore we simply return false.
    return false
end

end # module
