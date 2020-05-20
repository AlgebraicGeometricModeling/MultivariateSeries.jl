export ms_decompose, cst_rkf, eps_rkf, weights, normlz

import LinearAlgebra: diagm

#------------------------------------------------------------------------
eps_rkf = eps::Float64 -> function (S)
  i :: Int = 1;
  while i<= length(S) && S[i]/S[1] > eps
    i+= 1;
  end
  i-1;
end

cst_rkf = r::Int64 -> function (S) return r end

#------------------------------------------------------------------------
# Decomposition of the pencil of matrices
function ms_decompose(H::Vector{Matrix{C}}, lambda::Vector, rkf::Function) where C
    n = length(H)
    
    H0 = sum(H[i]*lambda[i] for i in 1:length(lambda))

    U, S, V = svd(H0)       # H0= U*diag(S)*V'
    r = rkf(S)

    Sr  = S[1:r]
    Sri = diagm([one(C)/S[i] for i in 1:r])

    M = Matrix{C}[]
    for i in 1:length(H)
    	push!(M, Sri*conj(U[:,1:r]')*H[i]*V[:,1:r])
    end

    Xi, E, DiagInfo = diagonalization(M)

    Uxi = (U[:,1:r].*Sr')*E
    Vxi = (E\ V[:,1:r]')

    Info = Dict{String,Any}( "diagonalization" => DiagInfo)
    return Xi, Uxi, Vxi, Info
end

#------------------------------------------------------------------------
"""
```
ms_decompose(σ :: Series{C,M}, rkf :: Function)
```
Decompose the series ``σ`` as a weighted sum of exponentials.
Return ``ω``, ``Ξ`` where
 - ``ω`` is the vector of weights,
 - ``Ξ`` is the matrix of frequency points, stored per row.
The list of monomials of degree ``\\leq {d-1 \\over 2}`` are used to construct
the Hankel matrix, where ``d`` is the maximal degree of the moments in ``σ``.

The optional argument `rkf` is the rank function used to determine the numerical rank from the vector S of singular values. Its default value `eps_rkf(1.e-6)` determines the rank as the first i s.t. S[i+1]/S[i]< 1.e-6 where S is the vector of singular values.

If the rank function cst_rkf(r) is used, the SVD is truncated at rank r.
"""
function ms_decompose(sigma::Series{R,M},
                      rkf::Function = eps_rkf(1.e-6),
                      weps::Float64=1.e-5) where {R, M}
    d = maxdegree(sigma)
    X = variables(sigma)

    d0 = div(d-1,2); d1 = d-1-d0
    B0 = monomials(X, seq(0:d0))
    B1 = monomials(X, seq(0:d1))

    H = Matrix{R}[hankel(sigma, B0, B1)]
    for x in X
        push!(H, hankel(sigma, B0, [b*x for b in B1]))
    end

    lambda = [1.0]
    Xi, Uxi, Vxi = ms_decompose(H, lambda,  rkf)

    n, r = size(Xi)
    
    w = fill(one(eltype(Xi)),r)

    for i in 1:r
        w[i] = Xi[1,i]
        Xi[:,i]/= Xi[1,i]
        w[i]*= Uxi[1,i]*Vxi[i,1]
    end

    # remove weights below threshold weps
    I = Bool[abs(w[i])>weps for i in 1:length(w)]
    return w[I], Xi[2:end,I]

end

#------------------------------------------------------------------------
function normlz(M::AbstractMatrix,i)
    diagm(0 => [1/M[i,j] for j in 1:size(M,1)])*M
end

