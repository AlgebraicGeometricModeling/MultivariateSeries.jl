
function norm_off(M)
    sqrt(norm(M)^2-norm(diag(M))^2)
end

function diagonalization_iter(D)
    n = size(D[1],1)
    s = length(D)
    
    X = fill(zero(D[1][1,1]),n,n)
    Y = fill(zero(D[1][1,1]),n,n)

    A = fill(zero(D[1][1,1]),s,2)
    b = fill(zero(D[1][1,1]),s)
    for i in 1:n
        for j in 1:n
            if i != j
                for k in 1:s
                    A[k,1] = D[k][i,i]
                    A[k,2] = D[k][j,j]
                    b[k]   = -D[k][i,j]
                end
                v = A\b
                X[i,j] =  v[1]
                Y[i,j] =  v[2]
            end
        end
    end
    for i in 1:n
        X[i,i]=1
        Y[i,i]=1
    end
    return X, Y
end

function diagonalization(M::Vector{Matrix{C}}, N=10, eps=1.e-3) where C
    n  = length(M)
    r  = size(M[1],1)

    M1 = sum(M[i]*randn(Float64) for i in 1:n)
    E  = eigvecs(M1)

    F  = inv(E)
    
    D  = vcat([Matrix{C}(I,r,r)],[F*M[i]*E for i in 1:length(M)])
    err = sum(norm_off.(D))
    delta = sum(norm.(D))
    println("diag off: ", err)

    if err/delta > 5.e-2
        nit = 0
        delta = err
        while nit < N && delta > eps
            err0 = err
            X,Y = diagonalization_iter(D)
            D = [Y*D[i]*X for i in 1:length(D)]
            E = E*X
            F = Y*F
            nit+=1
            err = sum(norm_off.(D))
            delta = err0-err
            #println("Off", nit,": ", err, "   delta: ", delta)
        end
        println("diag off: ", err, "  N: ",nit)
    end
    Xi = fill(zero(E[1,1]),n,r)
    for i in 1:r
    	for j in 1:n
	    Xi[j,i] = D[j+1][i,i]/D[1][i,i]
            #Xi[j,i] =(E[:,i]\(M[j]*E[:,i]))[1]
	end
    end
    return Xi, E
end

