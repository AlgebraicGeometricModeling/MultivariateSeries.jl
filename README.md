The package `MultivariateSeries.jl` provides tools for the manipulation of
series indexed by monomial exponents, sequence of moments, linear functionals on polynomials
and polynomial-exponential decomposition.

## Installation

To install the registered package within julia:

```julia
] add MultivariateSeries
```

or the latest version

```julia
] add https://github.com/AlgebraicGeometricModeling/MultivariateSeries.jl.git
```



## Example

```julia
using MultivariateSeries, DynamicPolynomials

X = @polyvar x1 x2 
n = length(X)
d = 4
r = 4

Xi0 = randn(n,r)
w0  = rand(r)

L = monomials(X,0:5)
sigma = series(w0, Xi0, L)


L2 = monomials(X,0:2)
L3 = monomials(X,0:3)
H = hankel(sigma, L2, L3)

w, Xi = decompose(sigma)
```

## Documentation
    
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://AlgebraicGeometricModeling.github.io/MultivariateSeries.jl/)

## More information

   - [source](https://github.com/AlgebraicGeometricModeling/MultivariateSeries.jl)
   - More information on [Julia](https://julialang.org/)


## Dependencies

- Julia 1.0
- DynamicPolynomials
- MultivariatePolynomials
