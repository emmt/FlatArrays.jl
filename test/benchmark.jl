module FlatArraysBenchmarking

using BenchmarkTools, FlatArrays
using Base:
    has_offset_axes

function print_left_justified(str::AbstractString, len::Integer, c::Char=' ')
    print(str)
    (n = len - length(str)) > 0 && print(c^n)
end

function print_right_justified(str::AbstractString, len::Integer, c::Char=' ')
    (n = len - length(str)) > 0 && print(c^n)
    print(str)
end

function runtests()
    dims = (3,4,5,6,7)
    A = randn(dims)
    B = view(A, :, :, :, 1:4, 3) # flat view
    C = view(A, :, :, :, 2:5, :) # non-flat view
    @assert isflat(A)
    @assert isflat(B)
    @assert !isflat(C)

    nelem = length(A)
    nrows = prod(dims[1:3])
    ncols = prod(dims[4:end])
    V = flatten(A, nelem)
    M = flatten(A, nrows, ncols)

    l = 30
    c = '.'

    println("In the following tests, A is a regular array, B and C are views of A,")
    println("B is flat, C is not flat, V is a flat vector from A and M is a flat")
    println("matrix from A.")
    println()
    print_left_justified("ndims(A) ",               l, c); @btime ndims($A);
    print_left_justified("ndims(B) ",               l, c); @btime ndims($B);
    print_left_justified("ndims(C) ",               l, c); @btime ndims($C);
    print_left_justified("ndims(V) ",               l, c); @btime ndims($V);
    print_left_justified("ndims(M) ",               l, c); @btime ndims($M);
    print_left_justified("firstindex(A) ",          l, c); @btime firstindex($A);
    print_left_justified("firstindex(B) ",          l, c); @btime firstindex($B);
    print_left_justified("firstindex(C) ",          l, c); @btime firstindex($C);
    print_left_justified("firstindex(V) ",          l, c); @btime firstindex($V);
    print_left_justified("firstindex(M) ",          l, c); @btime firstindex($M);
    print_left_justified("lastindex(A) ",           l, c); @btime lastindex($A);
    print_left_justified("lastindex(B) ",           l, c); @btime lastindex($B);
    print_left_justified("lastindex(C) ",           l, c); @btime lastindex($C);
    print_left_justified("lastindex(V) ",           l, c); @btime lastindex($V);
    print_left_justified("lastindex(M) ",           l, c); @btime lastindex($M);
    print_left_justified("iterate(A) ",             l, c); @btime iterate($A);
    print_left_justified("iterate(B) ",             l, c); @btime iterate($B);
    print_left_justified("iterate(C) ",             l, c); @btime iterate($C);
    print_left_justified("iterate(V) ",             l, c); @btime iterate($V);
    print_left_justified("iterate(M) ",             l, c); @btime iterate($M);
    print_left_justified("IndexStyle(A) ",          l, c); @btime IndexStyle($A);
    print_left_justified("IndexStyle(B) ",          l, c); @btime IndexStyle($B);
    print_left_justified("IndexStyle(C) ",          l, c); @btime IndexStyle($C);
    print_left_justified("IndexStyle(V) ",          l, c); @btime IndexStyle($V);
    print_left_justified("IndexStyle(M) ",          l, c); @btime IndexStyle($M);
    print_left_justified("has_offset_axes(A) ",     l, c); @btime has_offset_axes($A);
    print_left_justified("has_offset_axes(B) ",     l, c); @btime has_offset_axes($B);
    print_left_justified("has_offset_axes(C) ",     l, c); @btime has_offset_axes($C);
    print_left_justified("has_offset_axes(V) ",     l, c); @btime has_offset_axes($V);
    print_left_justified("has_offset_axes(M) ",     l, c); @btime has_offset_axes($M);

    print_left_justified("reshape(A,nelem) ",       l, c); @btime reshape($A,$nelem);
    print_left_justified("flatten(A,nelem) ",       l, c); @btime flatten($A,$nelem);
    print_left_justified("reshape(A,nrows,ncols) ", l, c); @btime reshape($A,$nrows,$ncols);
    print_left_justified("flatten(A,nrows,ncols) ", l, c); @btime flatten($A,$nrows,$ncols);
    print_left_justified("isflat(A) ",              l, c); @btime isflat($A);
    print_left_justified("isflat(B) ",              l, c); @btime isflat($B);
    print_left_justified("isflat(C) ",              l, c); @btime isflat($C);
    print_left_justified("isflat(V) ",              l, c); @btime isflat($V);
    print_left_justified("isflat(M) ",              l, c); @btime isflat($M);
    println()

end
runtests();
end # module
