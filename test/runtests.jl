#include("benchmark.jl")

module FlatArraysTests

using Test, FlatArrays, Random
using Base:
    IndexStyle,
    has_offset_axes

for d in 1:2,
    (f, g) in ((firstindex, :firstindex),
               (lastindex,  :lastindex),
               (size,       :size),
               (axes,       :axis),
               (stride,     :stride))
    @eval $(Symbol(g,d))(A) = $f(A, $d)
end

#function memcpy!(A::DenseArray{T}, i::Integer,
#                 B::DenseArray{T}, j::Integer, n::Integer) where {T}
#    @assert firstindex(A) ≤ i ≤ i + n - 1 ≤ lastindex(A)
#    @assert firstindex(B) ≤ j ≤ j + n - 1 ≤ lastindex(B)
#    ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
#          A, B, n*sizeof(T))
#end

function memcpy!(A::DenseArray{T}, B::DenseArray{T}, n::Integer) where {T}
    @assert firstindex(A) == 1 && 0 ≤ n ≤ lastindex(A)
    @assert firstindex(B) == 1 && 0 ≤ n ≤ lastindex(B)
    ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
          A, B, n*sizeof(T))
    return A
end

function samevalues(A::DenseArray{T}, B::DenseArray{T}) where {T}
    @assert length(A) == length(B)
    i = 0
    for a in A
        i += 1
        B[i] == a || return false
    end
    return true
end


@testset "FlatArrays" begin
    dims = (2,3,4,5)
    A = randn(dims)
    nelem = length(A)
    nrows = prod(dims[1:2])
    ncols = prod(dims[3:end])
    V = flatten(A, nelem)
    M = flatten(A, nrows, ncols)
    Q = ((firstindex,      1),
         (firstindex1,     1),
         (lastindex,       nelem),
         (lastindex1,      nelem),
         (length,          nelem),
         (ndims,           1),
         (eltype,          eltype(A)),
         (size,           (nelem,)),
         (size1,           nelem),
         (axes,           (1:nelem,)),
         (axis1,           1:nelem),
         (isflat,          true),
         (has_offset_axes, false))
    for i in randperm(length(Q))
        (f, r) = Q[i]
        @test f(V) == r
    end
    Q = ((firstindex,      1),
         (firstindex1,     1),
         (firstindex2,     1),
         (lastindex,       nelem),
         (lastindex1,      nrows),
         (lastindex2,      ncols),
         (length,          nelem),
         (ndims,           2),
         (eltype,          eltype(A)),
         (size,           (nrows,ncols)),
         (size1,           nrows),
         (size2,           ncols),
         (axes,           (1:nrows, 1:ncols)),
         (axis1,           1:nrows),
         (axis2,           1:ncols),
         (isflat,          true),
         (has_offset_axes, false))
    for i in randperm(length(Q))
        (f, r) = Q[i]
        @test f(M) == r
    end
    @test reshape(A, nelem) == V
    @test reshape(A, nrows, ncols) == M
    B = similar(A, nrows, ncols)
    @test memcpy!(B, M, nelem) == reshape(A, nrows, ncols)
    @test samevalues(B, M)
    @test samevalues(M, B)
    @test samevalues(B, V)
    @test samevalues(V, B)
end

end # module
