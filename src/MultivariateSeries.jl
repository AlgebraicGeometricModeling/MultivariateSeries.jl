module MultivariateSeries

using LinearAlgebra
using MultivariatePolynomials
using DynamicPolynomials

include("polynomials.jl")
include("series.jl")
include("moments.jl")
include("hankel.jl")
include("newton.jl")
include("decompose.jl")

export Seq, seq

mutable struct Seq{T} val::T end

function seq(args...)
    if length(args)>1
        Seq([args...])
    else
        Seq(args[1])
    end
end

end
