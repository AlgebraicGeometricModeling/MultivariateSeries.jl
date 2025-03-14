module MultivariateSeries

using LinearAlgebra
#using MultivariatePolynomials
using DynamicPolynomials

include("series.jl")
include("polynomials.jl")
include("moments.jl")
include("hankel.jl")
include("invsys.jl")

include("newton.jl")
include("diagonalisation.jl")
include("decompose.jl")

#=
export Seq, seq

mutable struct Seq{T} val::T end

function seq(args...)
    if length(args)>1
        Seq([args...])
    else
        Seq(args[1])
    end
end

function DynamicPolynomials.MonomialVector(V::Vector{DynamicPolynomials.Variable}, rg::Seq)
     L = DynamicPolynomials.Monomial{true}[]
     for i in rg.val
         append!(L, DynamicPolynomials.monomials(V,i))
     end
     L
end

function DynamicPolynomials.MonomialVector(V::Vector{DynamicPolynomials.Variable}, rg::UnitRange{Int64})
     L = DynamicPolynomials.Monomial{true}[]
     for i in rg
         append!(L, DynamicPolynomials.monomials(V,i))
     end
     L
end
=#

end
