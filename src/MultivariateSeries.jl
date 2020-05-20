module MultivariateSeries

using LinearAlgebra
using MultivariatePolynomials
using DynamicPolynomials

include("polynomials.jl")
include("series.jl")
include("moments.jl")
include("hankel.jl")
include("newton.jl")
include("diagonalisation.jl")
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

function DynamicPolynomials.MonomialVector(V::Vector{PolyVar{true}}, rg::Seq)
     L = DynamicPolynomials.Monomial{true}[]
     for i in rg.val
         append!(L, DynamicPolynomials.monomials(V,i))
     end
     L
end

end
