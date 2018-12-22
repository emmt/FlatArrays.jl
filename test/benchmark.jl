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

# Use eachindex()
function add0!(dest::AbstractArray{T,N},
               A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {T,N}
    # Checks must be skipped because they take time.
    #@assert axes(dest) == axes(A) == axes(B)
    @inbounds @simd for i in eachindex(dest, A, B)
        dest[i] = A[i] + B[i]
    end
    return dest
end

# Use linear indexing.
function add1!(dest::AbstractArray{T,N},
               A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {T,N}
    # Checks must be skipped because they take time.
    #@assert !Base.has_offset_axes(dest, A, B)
    #@assert (n = length(dest)) == length(A) == length(B)
    n = length(dest)
    @inbounds @simd for i in Base.OneTo(n)
        dest[i] = A[i] + B[i]
    end
    return dest
end

# Use 2D indexing.
function add2!(dest::AbstractArray{T,2},
               A::AbstractArray{T,2},
               B::AbstractArray{T,2}) where {T}
    # Checks must be skipped because they take time.
    #@assert (m = size(dest,1)) == size(A,1) == size(B,1)
    #@assert (n = size(dest,2)) == size(A,2) == size(B,2)
    m = size(dest,1)
    n = size(dest,2)
    @inbounds for j in 1:n
        @simd for i in 1:m
            dest[i,j] = A[i,j] + B[i,j]
        end
    end
    return dest
end

# Use 2D indexing + faking.
function add2p!(dest::AbstractArray{T,2},
                A::AbstractArray{T,2},
                B::AbstractArray{T,2}) where {T}
    # Checks must be skipped because they take time.
    #@assert (m = size(dest,1)) == size(A,1) == size(B,1)
    #@assert (n = size(dest,2)) == size(A,2) == size(B,2)
    m = size(dest,1)
    n = size(dest,2)
    @inbounds for j in 1:n
        off = m*(j - 1)
        @simd for i in 1:m
            k = i + off
            dest[k] = A[k] + B[k]
        end
    end
    return dest
end

# Use CartesianIndex()
function addI!(dest::AbstractArray{T,N},
               A::AbstractArray{T,N},
               B::AbstractArray{T,N}) where {T,N}
    # Checks must be skipped because they take time.
    #@assert axes(dest) == axes(A) == axes(B)
    @inbounds @simd for i in CartesianIndices(dest)
        dest[i] = A[i] + B[i]
    end
    return dest
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
    if false
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
    println("In the following tests, A is a regular array, R is a reshaped version of A,")
    println("M is a flat version of A and S is a random array of same size as R and M.")
    println()
    R = reshape(A, nrows, ncols)
    S = randn(size(R))
    print_left_justified("R'*R ",                   l, c); @btime $R'*$R;
    print_left_justified("S'*R ",                   l, c); @btime $S'*$R;
    print_left_justified("R'*S ",                   l, c); @btime $R'*$S;
    print_left_justified("M'*M ",                   l, c); @btime $M'*$M;
    print_left_justified("M'*R ",                   l, c); @btime $M'*$R;
    print_left_justified("R'*M ",                   l, c); @btime $R'*$M;
    println()
    end
    println("In the following tests:")
    println(" - `add0!` uses `eachindex`;")
    println(" - `add1!` uses 1D linear index;")
    println(" - `add2!` uses 2D linear indices;")
    println(" - `add2p!` uses 2D linear indices with a common stride for all arrays;")
    println(" - `addI!` uses CartesianIndex;")
    println()
    println("Tests with regular Array's of size $(dims):")
    A1 = randn(dims)
    A2 = randn(dims)
    A0 = similar(A1)
    print_left_justified("  add0!(A0, A1, A2) ", l, c); @btime add0!($A0, $A1, $A2);
    print_left_justified("  add1!(A0, A1, A2) ", l, c); @btime add1!($A0, $A1, $A2);
    print_left_justified("  addI!(A0, A1, A2) ", l, c); @btime addI!($A0, $A1, $A2);
    println()
    println("Tests with regular Vector's of length $(nelem):")
    B1 = randn(nelem) # true vector
    B2 = randn(nelem) # true vector
    B0 = similar(B1)  # true vector
    print_left_justified("  add0!(B0, B1, B2) ", l, c); @btime add0!($B0, $B1, $B2);
    print_left_justified("  add1!(B0, B1, B2) ", l, c); @btime add1!($B0, $B1, $B2);
    print_left_justified("  addI!(B0, B1, B2) ", l, c); @btime addI!($B0, $B1, $B2);
    println()
    println("Tests with FlatVector's of length $(nelem):")
    V1 = flatten(A1, nelem)
    V2 = flatten(A2, nelem)
    V0 = flatten(A0, nelem)
    print_left_justified("  add0!(V0, V1, V2) ", l, c); @btime add0!($V0, $V1, $V2);
    print_left_justified("  add1!(V0, V1, V2) ", l, c); @btime add1!($V0, $V1, $V2);
    print_left_justified("  addI!(V0, V1, V2) ", l, c); @btime addI!($V0, $V1, $V2);
    println()
    println("Tests with reshaped Vector's of length $(nelem):")
    R1 = reshape(A1, nelem)
    R2 = reshape(A2, nelem)
    R0 = reshape(A0, nelem)
    print_left_justified("  add0!(R0, R1, R2) ", l, c); @btime add0!($R0, $R1, $R2);
    print_left_justified("  add1!(R0, R1, R2) ", l, c); @btime add1!($R0, $R1, $R2);
    print_left_justified("  addI!(R0, R1, R2) ", l, c); @btime addI!($R0, $R1, $R2);
    println()
    println("Tests with regular Matrix's of size $((nrows, ncols)):")
    C1 = randn(nrows,ncols) # true matrix
    C2 = randn(nrows,ncols) # true matrix
    C0 = similar(C1)        # true matrix
    print_left_justified("  add0!(C0, C1, C2) ", l, c); @btime add0!($C0, $C1, $C2);
    print_left_justified("  add1!(C0, C1, C2) ", l, c); @btime add1!($C0, $C1, $C2);
    print_left_justified("  add2!(C0, C1, C2) ", l, c); @btime add2!($C0, $C1, $C2);
    print_left_justified("  add2p!(C0, C1, C2) ",l, c); @btime add2p!($C0, $C1, $C2);
    print_left_justified("  addI!(C0, C1, C2) ", l, c); @btime addI!($C0, $C1, $C2);
    println()
    println("Tests with FlatMatrix's of size $((nrows, ncols)):")
    M1 = flatten(A1, nrows, ncols)
    M2 = flatten(A2, nrows, ncols)
    M0 = flatten(A0, nrows, ncols)
    print_left_justified("  add0!(M0, M1, M2) ", l, c); @btime add0!($M0, $M1, $M2);
    print_left_justified("  add1!(M0, M1, M2) ", l, c); @btime add1!($M0, $M1, $M2);
    print_left_justified("  add2!(M0, M1, M2) ", l, c); @btime add2!($M0, $M1, $M2);
    print_left_justified("  add2p!(M0, M1, M2) ",l, c); @btime add2p!($M0, $M1, $M2);
    print_left_justified("  addI!(M0, M1, M2) ", l, c); @btime addI!($M0, $M1, $M2);
    println()
    println("Tests with reshaped Matrix's of size $((nrows, ncols)):")
    S1 = reshape(A1, nrows, ncols)
    S2 = reshape(A2, nrows, ncols)
    S0 = reshape(A0, nrows, ncols)
    print_left_justified("  add0!(S0, S1, S2) ", l, c); @btime add0!($S0, $S1, $S2);
    print_left_justified("  add1!(S0, S1, S2) ", l, c); @btime add1!($S0, $S1, $S2);
    print_left_justified("  add2!(S0, S1, S2) ", l, c); @btime add2!($S0, $S1, $S2);
    print_left_justified("  add2p!(S0, S1, S2) ",l, c); @btime add2p!($S0, $S1, $S2);
    print_left_justified("  addI!(S0, S1, S2) ", l, c); @btime addI!($S0, $S1, $S2);

end
runtests();
end # module
