#include("benchmark.jl")

module FlatArraysTests

using Test, FlatArrays, Random
using Base:
    IndexStyle,
    has_offset_axes

for d in 1:3,
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

function samevalues(A::AbstractArray{T}, B::AbstractArray{T},
                    exhaustive::Bool=false) where {T}
    @assert length(A) == length(B)
    j = 0
    for a in A
        j += 1
        B[j] == a || return false
    end
    if exhaustive
        i = 0
        for b in B
            i += 1
            A[i] == b || return false
        end
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
         (strides,        (1,)),
         (stride1,         1),
         (stride2,         nelem),
         (isflat,          true),
         (has_offset_axes, false),
         (IndexStyle,      IndexLinear()))
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
         (size,           (nrows, ncols)),
         (size1,           nrows),
         (size2,           ncols),
         (axes,           (1:nrows, 1:ncols)),
         (axis1,           1:nrows),
         (axis2,           1:ncols),
         (strides,        (1, nrows)),
         (stride1,         1),
         (stride2,         nrows),
         (stride3,         nelem),
         (isflat,          true),
         (has_offset_axes, false),
         (IndexStyle,      IndexLinear()))
    for i in randperm(length(Q))
        (f, r) = Q[i]
        @test f(M) == r
    end
    @test isflat(A) == true
    @test isflat(A) == true
    V1 = view(A, :,  :, 1:3, 3)  # flat view
    @test isflat(V1) == true
    V2 = view(A, :,  :, 2:3, :) # non-flat view
    @test isflat(V2) == false
    @test reshape(A, nelem) == V
    @test reshape(A, nrows, ncols) == M
    B = similar(A, nrows, ncols)
    @test memcpy!(B, M, nelem) == reshape(A, nrows, ncols)
    @test samevalues(B, M, true)
    @test samevalues(B, V, true)
    @test samevalues(flatten(V1,length(V1)), V1, true)
    @test samevalues(flatten(V2,length(V2)), V2, true)
end

end # module
