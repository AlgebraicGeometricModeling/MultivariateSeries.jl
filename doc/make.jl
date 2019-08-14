using Documenter
using LinearAlgebra
using MultivariatePolynomials
using MultivariateSeries

dir = "mrkd"
Expl = map(file -> joinpath("expl", file), filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"/code")))

makedocs(
    format = Documenter.HTML(),
    sitename = "MultivariateSeries",
    authors = "B. Mourrain",
    modules = [MultivariateSeries],
    build = "html",
    source = dir,
    pages = Any[
        "Home" => "index.md",
        "Example" => Expl,
        "Functions & types" => Code
    ],
    repo = "https://github.com/JuliaAlgebra/MultivariateSeries.jl/tree/master",
    doctest = false
)

deploydocs(
    repo = "github.com/JuliaAlgebra/MultivariateSeries.jl.git",
    target = "site",
    version = "v1.0.0",
    deps = nothing,
    make = nothing
)
