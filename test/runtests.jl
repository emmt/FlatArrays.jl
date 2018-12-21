module FlatArraysTests

using Test, FlatArrays, Random
#include("benchmark.jl")

@testset "FlatArrays" begin
    dims = (2,3,4,5)
    A = randn(dims)
    nelem = length(A)
    nrows = prod(dims[1:2])
    ncols = prod(dims[3:end])
    V = flatten(A, nelem)
    M = flatten(A, nrows, ncols)
    Q = ((firstindex, 1),
         (lastindex,  nelem),
         (length,     nelem),
         (ndims,      1),
         (eltype,     eltype(A)),
         (size,      (nelem,)),
         (axes,      (1:nelem,)))
    for i in randperm(length(Q))
        (f, r) = Q[i]
        @test f(V) == r
    end
    Q = ((firstindex, 1),
         (lastindex,  nelem),
         (length,     nelem),
         (ndims,      2),
         (eltype,     eltype(A)),
         (size,      (nrows,ncols)),
         (axes,      (1:nrows, 1:ncols)))
    for i in randperm(length(Q))
        (f, r) = Q[i]
        @test f(M) == r
    end
    @test reshape(A, nelem) == V
    @test reshape(A, nrows, ncols) == M
end

end # module
