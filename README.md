The package `MultivariateSeries.jl` provides tools for the manipulation of
series indexed by monomial exponents, sequence of moments, linear functionals on polynomials
and polynomial-exponential decomposition.

## Installation

To install the latest version of the package within julia:

```julia
] add https://github.com/AlgebraicGeometricModeling/MultivariateSeries.jl.git
```

or the registered package

```julia
] add MultivariateSeries
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

L = monomials(X,seq(0:5))
sigma = series(w0, Xi0, L)


L2 = monomials(X,0:2)
L3 = monomials(X,0:3)
H = hankel(sigma, L2, L3)

w, Xi = decompose(sigma)
```

## Documentation

   - [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://bmourrain.github.io/MultivariateSeries.jl/latest)
   - [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmourrain.github.io/MultivariateSeries.jl/dev)
   - More information on [Julia](https://julialang.org/)


## Dependencies

- Julia 1.0
- DynamicPolynomials
- MultivariatePolynomials
