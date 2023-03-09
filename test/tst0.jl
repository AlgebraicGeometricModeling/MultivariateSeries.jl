using MultivariateSeries

X = @ring x y z
n = length(X)
d = 4
r = 4

Xi = rand(n,r)
w = fill(1.0,r)
sigma = series(w,Xi,X, d)

k = 2
L = monomials(X,0:2)
H = hankel(sigma,L,L)

sigma1 = x^2*y*sigma
