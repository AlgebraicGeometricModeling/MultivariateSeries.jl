The package `MultivariateSeries.jl` provides tools for the manipulation of
series indexed by monomial exponents, sequence of moments, linear functionals on polynomials
and polynomial-exponential decomposition.

## Installation

To install the package within julia:

```julia
using Pkg
Pkg.clone("https://github.com/JuliaAlgebra/MultivariateSeries.jl.git")
```

## Example

```julia
using MultivariateSeries

X = @ring x1 x2 
n = length(X)
d = 4
r = 4

Xi0 = randn(n,r)
w0  = rand(r)

L = monoms(X,5)
sigma = series(w0, Xi0, L)


L2 = monoms(X,2)
L3 = monoms(X,3)
H = hankel(sigma, L2, L3)

w, Xi = decompose(sigma)
```

## Documentation

- [LATEST](http://www-sop.inria.fr/members/Bernard.Mourrain/software/MultivariateSeries/index.html)
- More information on [Julia](https://julialang.org/)


## Dependencies

- Julia 1.0
- DynamicPolynomials
- MultivariatePolynomials
