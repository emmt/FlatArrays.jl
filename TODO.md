* Extend to more than 1 or 2 dimensional arrays.  Use `@generated` for number
  of dimensions > 2.

* `flatten` should have a *lazy* counterpart which just returns its
  argument if already suitable.

* `flatten(A)` with no dimensions yields a *flat* array of same dimensions
  as `A`.

* Implement `checkbound` and other related methods on FlatArray's.

* Allow for a mixture of `CartesianIndex` and integers.
