###
### Series as linear functionals on polynomials.
###
### Hankel matrices associated to series.
###
### Bernard Mourrain
###

export series, hankel, hankelbasis, matrixof

#-----------------------------------------------------------------------
"""
```
hankel(σ::Series{C,M}, L1::Vector{M}, L2::Vector{M}) -> Array{C,2}
```
Hankel matrix of ``σ`` with the rows indexed by the list of
polynomials L1 and the columns by L2. The entries are the dot product for ``σ``
of the corresponding elements in L1 and L2.

Example
-------
```
julia> L =[1, x1, x2, x1^2, x1*x2, x2^2]

julia> H = hankel(s,L,L)
6x6 Array{Float64,2}:
  4.0   5.0   7.0    5.0  11.0  13.0
  5.0   5.0  11.0   -1.0  17.0  23.0
  7.0  11.0  13.0   17.0  23.0  25.0
  5.0  -1.0  17.0  -31.0  23.0  41.0
 11.0  17.0  23.0   23.0  41.0  47.0
 13.0  23.0  25.0   41.0  47.0  49.0
```
"""
function hankel(sigma::Series{C,M}, L1::AbstractVector, L2::AbstractVector) where {C, M<:AbstractMonomial}
    return [sigma[L1[i]*L2[j]] for i in 1:length(L1), j in 1:length(L2)]
end

#------------------------------------------------------------------------
"""
```
series(H::Matrix{C}, L1::Vector{M}, L2::Vector{M}) -> Series{C,M}
```
Compute the series associated to the Hankel matrix H, with rows (resp. columns)
indexed by the array of monomials L1 (resp. L2).
"""
function series(H::Array{C,2} , L1::Vector{M}, L2::Vector{M}) where {C, M<:AbstractMonomial}
    res = zero(Series{C,M})
    i = 1
    for m1  in L1
        j=1
        for m2 in L2
            m = m1*m2
            if res[m] == zero(C)
	        res[m] = H[i,j]
	    end
            j+=1
        end
        i+=1
    end
   res
end

#----------------------------------------------------------------------
"""
Generate the table of Hankel matrices ``H_α`` associated to the monomials ``x^α`` generating the space of Hankel matrices indexed by L1 (for the rows) and L2 (for the columns).
"""
function hankelbasis(L1,L2)
    m1 = length(L1)
    m2 = length(L2)
    table = Dict( L1[1]*L2[1] => fill(0.,m1,m2))
    for i in 1:m1
        for j in 1:m2
            M = get(table,L1[i]*L2[j], 0)
            if M != 0
                M[i,j]=1.0
            else
                M = fill(0.,m1, m2)
                M[i,j]=1.0
                table[L1[i]*L2[j]] = M
            end
        end
    end
    table
end


"""
```
matrixof(P::Vector,L::Vector)
```
  matrix `M` which rows are the coefficient the polynomials in P with respect to the monomial vector L, so that `M*L=P` 
"""
function matrixof(P::AbstractVector, L::AbstractVector) 

    M = fill(zero(coefficient_type(P[1])), length(P), length(L))
    idx = Dict{typeof(L[1]),Int64}()
    for i in 1:length(L)
        idx[L[i]] = i
    end

    for i in 1:length(P)
        for (c,m) in zip(coefficients(P[i]),monomials(P[i]))
            j = get(idx,m,0)
            if j!=0 M[i,j]=c end
        end
    end
    M
end


"""
```
M, L = matrixof(S::Vector{Series})
```
 -  `M` matrix which rows are the coefficient the polynomials in P with respect to the monomial vector L
 -  `L` vector of all monomials of `S`
"""

function matrixof(S::Vector{Series{C,M}}) where {C,M}

    L = monomials(sum(p for p in [sum(monomials(s)) for s in S]))
    mtr = fill(zero(C), length(S), length(L) )
    for i in 1:length(L)
        for j in 1: length(S)
            mtr[j,i] = S[j][L[i]]
        end
    end
    return mtr, L
end




    
