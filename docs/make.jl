using Documenter
using LinearAlgebra
using MultivariatePolynomials
using MultivariateSeries

dir = joinpath(pwd(),"mrkd")
Expl = map(file -> joinpath("expl", file), filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"/code")))

makedocs(
    sitename = "MultivariateSeries",
    format = Documenter.HTML(),
    authors = "B. Mourrain",
    modules = [MultivariateSeries],
    build = "MultivariateSeries.jl/docs",
    source = dir,
    pages = Any[
        "Home" => "index.md",
        "Example" => Expl,
        "Functions & types" => Code
    ],
    doctest = false
)

deploydocs(
    repo = "github.com/AlgebraicgeometricModeling/MultivariateSeries.jl.git"
)
