# FlatArrays

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/FlatArrays.jl.svg?branch=master)](https://travis-ci.org/emmt/FlatArrays.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/FlatArrays.jl?branch=master)](https://ci.appveyor.com/project/emmt/FlatArrays-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/emmt/FlatArrays.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/emmt/FlatArrays.jl?branch=master)
[![codecov.io](http://codecov.io/github/emmt/FlatArrays.jl/coverage.svg?branch=master)](http://codecov.io/github/emmt/FlatArrays.jl?branch=master)

This [Julia](http://julialang.org/) package provides *flat* arrays.  That is
arrays with 1-based indices, whose elements are all continuous and in
column-major order.  Such arrays are useful for writting efficient Julia code
using linear indexing.  This kind of arrays is also what is expected by most
numerical libraries and FlatArrays may be helpful to link such libraries in
Julia.

FlatArrays are built from any Julia array and can have any dimensions providing
that their total size do not exceed that of the input array.  This is like
`reshape` except that building a `FlatArray` from a Julia `Array` is faster
than calling `reshape` (by a factor ~ 7 for a vector, ~ 5 for a matrix).  This
may be significant for very small arrays.

Indexing of FlatArray's is as fast as regular Julia arrays.

I found FlatArrays useful to pretend that multidimensional arrays are vectors
or matrices so as to automatically benefit from very fast linear algebra
methods (*e.g.*, BLAS or LAPACK).  Again, this seems like `reshape` except that
`reshape` is slower and does not guarantee how the elements are stored.


## Simple Usage

Import the package:

```julia
using FlatArrays
```

and then:

```julia
flatten(A, n) -> V
```

yields a *flat* vector `V` of length `n` whose elements are given by array `A`.
To build a *flat* matrix from the elements of `A`, call:

```julia
flatten(A, m, n) -> M
```

which yields a *flat* matrix `M` of size `m` rows by `n` columns whose elements are
given by array `A`.


## Installation

To install FlatArrays, you must clone its repository:

```sh
using Pkg
Pkg.clone("https://github.com/emmt/FlatArrays.jl.git")
```

The FlatArrays package is pure Julia code and there is nothing to build.
