using MultivariateSeries
X = @ring x1 x2
d = 4

w = [1, 1.5, -2.0]
Xi = [ 1 1;
       0 0;
       0 -1.0]'

F = series(w, Xi, X, d)

decompose(F)


