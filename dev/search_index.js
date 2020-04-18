var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#MultivariateSeries-1",
    "page": "Home",
    "title": "MultivariateSeries",
    "category": "section",
    "text": "Package for the decomposition of tensors and polynomial-exponential series."
},

{
    "location": "#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "The package MultivariateSeries.jl provides tools for the manipulation of sequences (sigma_alpha)_alpha in mathbbK^mathbbN^n indexed by multivariate indices alpha in mathbbN^n which are represented as series:      sigma(mathbfz) = sum_alpha in mathbbN^n sigma_alpha mathbfz^alphaSometimes it is preferrable to associate the following series in mathbfy to the sequence sigma:sigma(mathbfy) = sum_alpha in mathbbN^n sigma_alpha fracmathbfy^alphaalphaThe sequence sigma or the series sigma(z) represents a linear functional on the polynomials:     sigma p= sum_alpha in mathbbN^n p_alpha mathbfx^alpha mapsto sum_alpha in mathbbN^n p_alpha sigma_alphaThe series are represented as association tables (or dictionnary) between (a finite set of) monomials and coefficients. They are printed using dual variables dxi:using MultivariateSeries\nX = @ring x1 x2\n\njulia> p = (1+x1)^3 + 0.5*(1+x1+x2)^2\nx1³ + 3.5x1² + x1x2 + 0.5x2² + 4.0x1 + x2 + 1.5\n\njulia> sigma = dual(p)\ndx1^3 + 1.5 + 3.5dx1^2 + 4.0dx1 + dx2 + dx1*dx2 + 0.5dx2^2Since the series is represented by a table, the order in which the dual monomials are printed is the order used in the table. It is not necessarily sorted by a monomial ordering.julia> sigma.terms\nDict{Monomial{true},Float64} with 7 entries:\n   x1³  => 1.0\n   1    => 1.5\n   x1²  => 3.5\n   x1   => 4.0\n   x2   => 1.0\n   x1x2 => 1.0\n   x2²  => 0.5Series act as linear functionals on polynomials via the dot product:julia> dot(sigma,x1^2)\n3.5\njulia> dot(sigma,x2^4)\n0.0"
},

{
    "location": "#Polynomial-exponential-decomposition-1",
    "page": "Home",
    "title": "Polynomial-exponential decomposition",
    "category": "section",
    "text": "The package provide tools for solving the following decomposition problem:Given (the first terms of) sequence sigma in mathbbK^mathbbN^n or the series  sigma(mathbfy) in mathbbKmathbfy, we want to decompose it as polynomial-exponential series sigma(mathbfy) = sum_i=1^r omega_i(mathbfy) e^xi_i1 y_1+ cdots + xi_in y_nwith polynomials omega_i(mathbfy) and points xi_i= (xi_i1 ldots xi_in)in mathbbK^n.  omega_i are called the weights and  xi_i the frequencies of the decomposition.These types of decompositions appear in many problems (see Examples). The package MultivariateSeries provides functions to manipulate (truncated) series, to construct truncated Hankel matrices, and to compute such a decomposition from these Hankel matrices."
},

{
    "location": "#sec_examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "Pages = map(file -> joinpath(\"expl\", file), filter(x ->endswith(x, \"md\"), readdir(\"expl\")))"
},

{
    "location": "#Functions-and-types-1",
    "page": "Home",
    "title": "Functions and types",
    "category": "section",
    "text": "Pages = map(file -> joinpath(\"code\", file), filter(x ->endswith(x, \"md\"), readdir(\"code\"))) "
},

{
    "location": "#sec_installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is available at https://github.com/bmourrain/MultivariateSeries.jl.gitTo install it from Julia:using Pkg\nPkg.clone(\"https://github.com/bmourrain/MultivariateSeries.jl.git\")It can then be used as follows:using MultivariateSeriesSee the Examples for more details."
},

{
    "location": "#Dependencies-1",
    "page": "Home",
    "title": "Dependencies",
    "category": "section",
    "text": "The package MultivariateSeries depends on the following packages:DynamicPolynomials package on multivariate polynomials represented as lists of monomials.\nMultivariatePolynomials generic interface package for multivariate polynomials.These packages will be installed with MultivariateSeries  (see installation)."
},

{
    "location": "expl/0.Decomposition/#",
    "page": "Decomposition algorithm",
    "title": "Decomposition algorithm",
    "category": "page",
    "text": ""
},

{
    "location": "expl/0.Decomposition/#Decomposition-algorithm-1",
    "page": "Decomposition algorithm",
    "title": "Decomposition algorithm",
    "category": "section",
    "text": "using MultivariateSeries\nusing LinearAlgebra\nX = @ring x1 x2 x33-element Array{DynamicPolynomials.PolyVar{true},1}:\n x1\n x2\n x3We want to find a sparse representation of the following series known up to degree 3:sigma = dual(6.0 + 4.0*x1 + 15.0*x2 + 6.0*x3 + 6.0*x1^2 + 20.0*x1*x2 + 4.0*x1*x3 + 43.0*x2^2 + 15.0*x2*x3 + 6.0*x3^2 - 26.0*x1^3 + 30.0*x1^2*x2 + 6.0*x1^2*x3 + 72.0*x1*x2^2 + 20.0*x1*x2*x3 + 4.0*x1*x3^2 + 129.0*x2^3 + 43.0*x2^2*x3 + 15.0*x2*x3^2 + 6.0*x3^3)4.0dx1*dx3 + 15.0dx2 + 6.0dx3 + 20.0dx1*dx2 + 6.0dx3^3 + 43.0dx2^2dx3 - 26.0dx1^3 + 129.0dx2^3 + 30.0dx1^2dx2 + 15.0dx2*dx3 + 20.0dx1*dx2*dx3 + 6.0dx1^2 + 6.0dx3^2 + 4.0dx1 + 43.0dx2^2 + 6.0dx1^2dx3 + 4.0dx1*dx3^2 + 72.0dx1*dx2^2 + 15.0dx2*dx3^2 + 6.0L1 = monoms(X,1)\nL2 = monoms(X,2)10-element Array{DynamicPolynomials.Monomial{true},1}:\n 1   \n x1  \n x2  \n x3  \n x1² \n x1x2\n x1x3\n x2² \n x2x3\n x3²H = hankel(sigma,L1,L2)4×10 Array{Float64,2}:\n  6.0   4.0  15.0   6.0    6.0  20.0   4.0   43.0  15.0   6.0\n  4.0   6.0  20.0   4.0  -26.0  30.0   6.0   72.0  20.0   4.0\n 15.0  20.0  43.0  15.0   30.0  72.0  20.0  129.0  43.0  15.0\n  6.0   4.0  15.0   6.0    6.0  20.0   4.0   43.0  15.0   6.0The rank of H_sigma will give us an idea on the dimension of mathcalA_sigma.rank(H)3We check that 1 x_1 x_2 is a basis of mathcalA_sigma: B0 = L1[1:3]3-element Array{DynamicPolynomials.Monomial{true},1}:\n 1 \n x1\n x2H0 = hankel(sigma, B0, B0)3×3 Array{Float64,2}:\n  6.0   4.0  15.0\n  4.0   6.0  20.0\n 15.0  20.0  43.0rank(H0)3Let us compute the shifted (truncated) Hankel operators.H1 = hankel(sigma, B0, B0*x1)\nH2 = hankel(sigma, B0, B0*x2)\nH3 = hankel(sigma, B0, B0*x3);\nH  = [H1,H2,H3]\nH[1]3×3 Array{Float64,2}:\n  4.0    6.0  20.0\n  6.0  -26.0  30.0\n 20.0   30.0  72.0M = [ H0^(-1)*H[i] for i in 1:3 ]\nM[1]3×3 Array{Float64,2}:\n  1.11022e-16   9.14286  -0.571429\n  1.0           3.85714   1.57143 \n -1.11022e-16  -4.28571   1.14286The eigenvalues and eigenvectors of M_x_1 areWe deduce the operators of multiplication by the variables in the basis B_0:v, E = eigen(M[1])Eigen{Float64,Float64,Array{Float64,2},Array{Float64,1}}\neigenvalues:\n3-element Array{Float64,1}:\n -0.9999999999999991\n  4.000000000000002 \n  2.000000000000002 \neigenvectors:\n3×3 Array{Float64,2}:\n  0.963087  -0.811107  -0.762001\n -0.120386  -0.324443  -0.127   \n -0.240772   0.486664   0.635001The matrices M_x_i are diagonal in this basis:D = [E^(-1)*M[i]*E for i in 1:3]\nD[1]3×3 Array{Float64,2}:\n -1.0          -6.99441e-15  -3.66374e-15\n  4.21885e-15   4.0          -4.44089e-15\n -4.66294e-15  -3.9968e-15    2.0D[2]3×3 Array{Float64,2}:\n  1.0          -4.44089e-16  -1.44329e-15\n  8.88178e-16   2.0           2.66454e-15\n -3.55271e-15   2.66454e-15   3.0D[3]3×3 Array{Float64,2}:\n  1.0           3.33067e-16  1.11022e-16\n -9.4369e-16    1.0          6.66134e-16\n  5.55112e-16  -6.66134e-16  1.0Looking at the corresponding terms on the diagonal, we get the coordinates of the points Xi:Xi = [ D[i][j,j] for i in 1:3, j in 1:3]3×3 Array{Float64,2}:\n -1.0  4.0  2.0\n  1.0  2.0  3.0\n  1.0  1.0  1.0We normalize the eigenvectors by v_i over v_i(xi_i) and get the interpolation polynomials at the points xi_i:Dg = E\'*vcat(fill(1.,1,3), Xi[1:2,:])\nE = E*Dg^(-1)\nU = E\'*B03-element Array{DynamicPolynomials.Polynomial{true,Float64},1}:\n -0.14285714285714324x1 - 0.2857142857142862x2 + 1.142857142857143 \n 0.28571428571428614x1 - 0.4285714285714279x2 + 0.7142857142857121 \n -0.14285714285714332x1 + 0.7142857142857134x2 - 0.8571428571428543We deduce the weights w_i=sigma(u_i):w = hankel(sigma, U, [L1[1]])3×1 Array{Float64,2}:\n  1.999999999999992 \n -1.0000000000000018\n  5.000000000000002Using the command decompose, we can get directly the same decomposition: w, Xi = decompose(sigma)([-1.0, 5.0, 2.0], [4.0 2.0 -1.0; 2.0 3.0 1.0; 1.0 1.0 1.0])Xi3×3 Array{Float64,2}:\n 4.0  2.0  -1.0\n 2.0  3.0   1.0\n 1.0  1.0   1.0w3-element Array{Float64,1}:\n -1.0000000000000129\n  5.000000000000011 \n  1.9999999999999998The series decomposes as 2 mathfrake_(-111) + 5 mathfrake_(231) - mathfrake_(421)."
},

{
    "location": "expl/3.MultivariateProny/#",
    "page": "Multivariate exponential decompositon",
    "title": "Multivariate exponential decompositon",
    "category": "page",
    "text": ""
},

{
    "location": "expl/3.MultivariateProny/#Multivariate-exponential-decompositon-1",
    "page": "Multivariate exponential decompositon",
    "title": "Multivariate exponential decompositon",
    "category": "section",
    "text": "using MultivariateSeriesWe consider the following function, which is a sum of 6 complex exponentials f = (u,v) -> 0.5*cos(0.7*pi*(u+v))+0.6*sin(4*pi*u)-0.2*cos(pi*v);(Image: waves)In order to recover the frequencies or exponents of these exponential terms and their coefficients, we sample the function on a grid (alpha_1 over T alpha_2 over T) alpha=(alpha_1alpha_2)in Asubset mathbbN^2. This defines a sequence of moments sigma_alpha=f(alpha_1 over T alpha_2 over T). We compute its generating series truncated in degree leq 5.X = @ring x1 x2\nL = monoms(X,5)\nT = 10\nmnt = (V->f(V[1]/T,V[2]/T))\nsigma = series(mnt, L)0.2569085959993554dx2^4 - 0.22417045976016967dx1^3dx2 + 0.7358257607718759dx1*dx2^4 + 0.8585922907464658dx1 + 0.5575373543042985dx1^2dx2 + 0.46210935078676274dx1^2dx2^3 + 0.6050846776084937dx1^2 + 0.02699524986977328dx1^5 - 0.15759364518763863dx1^3 + 0.2775204557293506dx2^3 + 0.29774707771034303dx2 - 0.2874793003806999dx1^3dx2^2 + 0.5095797473748392dx1^2dx2^2 + 0.2906101273580203dx2^2 + 0.7717888541929423dx1*dx2^3 + 0.22699524986977343dx2^5 - 0.5338499631663495dx1^4dx2 + 0.8328361327510712dx1*dx2 + 0.8039080170899477dx1*dx2^2 - 0.4519219149027473dx1^4 + 0.3Computing its decomposition using svdw, Xi = decompose(sigma);yields the weights omega of the exponential terms in f and the exponentials Xi:log.(Xi\')*T/pi6×2 Array{Complex{Float64},2}:\n 6.84243e-12+1.89601e-12im  -1.46851e-12-1.0im        \n 6.84243e-12-1.89601e-12im  -1.46851e-12+1.0im        \n 1.88616e-12-0.7im           4.84085e-13-0.7im        \n 1.88616e-12+0.7im           4.84085e-13+0.7im        \n 8.94708e-13-4.0im          -2.79182e-13+8.46664e-15im\n 8.94708e-13+4.0im          -2.79182e-13-8.46664e-15imBy taking the log and scaling by Tover pi, we recover the frequency vectors within precision 1O^-11. w6-element Array{Complex{Float64},1}:\n  0.08650714237109196 + 0.05016487136224489im\n  0.08650714237109196 - 0.05016487136224486im\n -0.17203277153436108 - 0.18139659731722088im\n  -0.1720327715343611 + 0.18139659731722083im\n -0.16999220553356456 - 0.24718950232123527im\n -0.16999220553356453 + 0.2471895023212353im"
},

{
    "location": "expl/4.DiracMeasure/#",
    "page": "Weighted sum of Dirac Measures",
    "title": "Weighted sum of Dirac Measures",
    "category": "page",
    "text": ""
},

{
    "location": "expl/4.DiracMeasure/#Weighted-sum-of-Dirac-Measures-1",
    "page": "Weighted sum of Dirac Measures",
    "title": "Weighted sum of Dirac Measures",
    "category": "section",
    "text": "using MultivariateSeriesSeries with 3 variablesx = @ring x1 x2 x3\nn = length(x)\nr = 4;Random weights in 01w0 = rand(Float64,r)4-element Array{Float64,1}:\n 0.10310570990285539\n 0.659237146402671  \n 0.42854695858385483\n 0.9028954619877378Random points in 01^nXi0 = rand(Float64,n,r)3×4 Array{Float64,2}:\n 0.894353  0.723909  0.847517   0.243868\n 0.621765  0.474676  0.889351   0.234888\n 0.448695  0.324577  0.0910113  0.929199Moment function of the sum of the Dirac measures of the points Xi_0 with weights omega_0 and its generating series up to degree 3.mt = moment(w0, Xi0)\ns = series(mt, monoms(x, 3))0.4339257411626396dx1*dx3 + 0.9702397337601841dx2 + 1.1382083719944895dx3 + 0.6585957641575078dx1*dx2 + 0.756553791183631dx3^3 + 0.14323371231850154dx2^2dx3 + 0.5978228622857406dx1^3 + 0.40844393047707755dx2^3 + 0.36208356970237654dx2*dx3 + 0.5016353149729792dx1^2dx2 + 0.17670723666808152dx1*dx2*dx3 + 0.7894555317442516dx1^2 + 0.8733277036829585dx3^2 + 1.1528284153895243dx1 + 0.5771697834203546dx2^2 + 0.22704558514092504dx1^2dx3 + 0.2619613422143817dx1*dx3^2 + 0.4425968864300194dx1*dx2^2 + 0.2321415223320108dx2*dx3^2 + 2.093785276877119Decomposition of the series from its terms up to degree 3.w, Xi = decompose(s);w4-element Array{Float64,1}:\n 0.4285469585838939 \n 0.10310570990277962\n 0.6592371464026707 \n 0.9028954619877747Xi3×4 Array{Float64,2}:\n 0.847517   0.894353  0.723909  0.243868\n 0.889351   0.621765  0.474676  0.234888\n 0.0910113  0.448695  0.324577  0.929199"
},

{
    "location": "expl/5.SparseInterpolation/#",
    "page": "Sparse interpolation",
    "title": "Sparse interpolation",
    "category": "page",
    "text": ""
},

{
    "location": "expl/5.SparseInterpolation/#Sparse-interpolation-1",
    "page": "Sparse interpolation",
    "title": "Sparse interpolation",
    "category": "section",
    "text": "using MultivariateSeriesA sparse polynomial in 3 variablesX = @ring x1 x2 x3\nf = 6.7x1^4*x2^5*x3 + 10.2x1^2*x2*x3^3 - 3.4x1*x2^2*x36.7x1^{4}x2^{5}x3 + 10.2x1^{2}x2x3^{3} - 3.4x1x2^{2}x3The series of moments f(zeta^alpha) for alphaleq 3.zeta = fill(0.9, length(X))\nsigma = series(f, zeta, X,3)7.2252810000000025dx1*dx3 + 10.382283000000001dx2 + 10.405800000000001dx3 + 7.5529172763000005dx1*dx2 + 6.357388987800002dx3^3 + 6.117862993803001dx2^2dx3 + 4.834376094422702dx1^3 + 7.008371185034148dx2^3 + 7.774274700000001dx2*dx3 + 5.495308104980432dx1^2dx2 + 5.5261037486700015dx1*dx2*dx3 + 6.8223503070000024dx1^2 + 8.0936982dx3^2 + 9.59787dx1 + 8.36740554867dx2^2 + 4.995745656300002dx1^2dx3 + 5.472820242000001dx1*dx3^2 + 6.217299094482389dx1*dx2^2 + 5.852477610000002dx2*dx3^2 + 13.499999999999998Computing its decomposition using svdw, Xi = decompose(sigma);yields the coefficients of the terms of f as the weights omega, and the exponents of the monomials of f as the log_zeta of the points Xi:w3-element Array{Float64,1}:\n -3.399999999998614\n 10.199999999997454\n  6.700000000001163Ex = log(Xi, zeta)3×3 Array{Float64,2}:\n 1.0  2.0  4.0\n 2.0  1.0  5.0\n 1.0  3.0  1.0"
},

{
    "location": "expl/6.Decoding/#",
    "page": "Decoding algebraic codes (BMS)",
    "title": "Decoding algebraic codes (BMS)",
    "category": "page",
    "text": ""
},

{
    "location": "expl/6.Decoding/#Decoding-algebraic-codes-(BMS)-1",
    "page": "Decoding algebraic codes (BMS)",
    "title": "Decoding algebraic codes (BMS)",
    "category": "section",
    "text": "We consider the code C formed by the words min mathbbk^l such that  $ m.[f(P1), \\ldots, f(Pl)]=0 $  where fin Vsubset mathbbkx_1x_n.using MultivariateSeries\nX = @ring x1 x22-element Array{DynamicPolynomials.PolyVar{true},1}:\n x1\n x2We consider the following points P and the vector space V spanned by the monomials M of degree le 2 in two variables x_1 x_2.P = [\n1   1;\n1  -1;\n-1  1;\n-1 -1;\n0   1; \n2  -1;\n1   2;\n1  -2]\'\n\nM = monoms(X,2)6-element Array{DynamicPolynomials.Monomial{true},1}:\n 1   \n x1  \n x2  \n x1² \n x1x2\n x2²The words of the code are the kernel of the following matrix:using DynamicPolynomials\nfunction vdm(P,L)\n    [ Polynomial{true,Float64}(L[j])(P[:,i])\n        for  j in 1:length(L),i in 1:size(P,2) ]\nend\nW = vdm(P,M)6×8 Array{Float64,2}:\n 1.0   1.0   1.0   1.0  1.0   1.0  1.0   1.0\n 1.0   1.0  -1.0  -1.0  0.0   2.0  1.0   1.0\n 1.0  -1.0   1.0  -1.0  1.0  -1.0  2.0  -2.0\n 1.0   1.0   1.0   1.0  0.0   4.0  1.0   1.0\n 1.0  -1.0  -1.0   1.0  0.0  -2.0  2.0  -2.0\n 1.0   1.0   1.0   1.0  1.0   1.0  4.0   4.0We receive the following word:r = [3, 3, 3, 0, -6, -2, 0, -1]8-element Array{Int64,1}:\n  3\n  3\n  3\n  0\n -6\n -2\n  0\n -1It is not a word of code C, since the following vector of syndroms is not zero:s = W*r6-element Array{Float64,1}:\n  0.0\n -2.0\n  1.0\n  0.0\n  3.0\n -3.0We want to correct it. For that, we build the corresponding series of syndroms:sigma = dual(M\'*s)3.0dx1*dx2 - 3.0dx2^2 + dx2 - 2.0dx1The Hankel matrix of sigma in degree le 1 is:L1 = monoms(X,1)\nH = hankel(sigma, L1, L1)3×3 Array{Float64,2}:\n  0.0  -2.0   1.0\n -2.0   0.0   3.0\n  1.0   3.0  -3.0An element in its kernel gives an error locator polynomial of degree 1:using LinearAlgebra\nle = nullspace(H); le/=le[3]\nple = (L1\'*le)[1]0.5000000000000002x1 + x2 + 1.5We check for which point in P, this polynomial vanishes. This will give the position where an error occurs:er = le\'*vcat(fill(1.,1,size(P,2)),P)1×8 Array{Float64,2}:\n 3.0  1.0  2.0  -2.22045e-16  2.5  1.5  4.0  0.0ie = []\nfor i in 1:length(er)\n    if isapprox(er[i],0.0;atol=1e-10) push!(ie, i) end\nend\nie2-element Array{Any,1}:\n 4\n 8These are the following points of P:E = ([P[j,ie[i]] for i in 1:length(ie), j in 1:size(P,1)])\'2×2 Adjoint{Int64,Array{Int64,2}}:\n -1   1\n -1  -2To get the error, that is the weights, we solve the system: E*omega =sigma_x_1 sigma_x_2:cr = E\\(W*r)[2:3]2-element Array{Float64,1}:\n  1.0\n -1.0We can now correct the received message by removing the weights cr at the positions of the errors ie:c=copy(r)\nfor i in 1:length(ie) c[ie[i]]-= cr[i] end \nc8-element Array{Int64,1}:\n  3\n  3\n  3\n -1\n -6\n -2\n  0\n  0We check that the corrected message is a word of the code:W*c6-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0"
},

{
    "location": "code/1.series/#",
    "page": "Series",
    "title": "Series",
    "category": "page",
    "text": ""
},

{
    "location": "code/1.series/#MultivariateSeries.Series",
    "page": "Series",
    "title": "MultivariateSeries.Series",
    "category": "type",
    "text": "Series{C,M}\n\nClass representing multivariate series. The series is a dictionary, which associates values of type C to monomials of type M.\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.dual",
    "page": "Series",
    "title": "MultivariateSeries.dual",
    "category": "function",
    "text": "dual(p::Polynomial) -> Series{C,M}\n\nCompute the series associated to the polynomial p, replacing the variables xi by their dual variables dxi. C is the type of coefficients  of the polynomial p and M its type of monomials.\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#LinearAlgebra.dot",
    "page": "Series",
    "title": "LinearAlgebra.dot",
    "category": "function",
    "text": "dot(σ::Series{C,M}, p::Monomial) -> C\ndot(σ::Series{C,M}, p::Term) -> C\ndot(σ::Series{C,M}, p::Polynomial) -> C\ndot(σ::Series{C,M}, p::Polynomial, q::Polynomial) -> C\n\nCompute the dot product  p q _σ =  σ  p q  or   σ  p  for p, q polynomials, terms or monomials. Apply the linear functional sigma on monomials, terms, polynomials \n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.moment",
    "page": "Series",
    "title": "MultivariateSeries.moment",
    "category": "function",
    "text": "moment(w::Vector{C}, P::Matrix{C}) -> Vector{Int64} -> C\n\nCompute the moment function α - _i ω_i P_i^α associated to the sequence P of r points of dimension n, which is a matrix of size r*n and the weights w.\n\n\n\n\n\nmoment(p::Polynomial, zeta::Vector{C}) -> Vector{Int64} -> C\n\nCompute the moment function α rightarrow p(ζ^α).\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.series",
    "page": "Series",
    "title": "MultivariateSeries.series",
    "category": "function",
    "text": "Construct the series with the term (c,m).\n\n\n\n\n\nConstruct the series from an array of pairs  m=>c where m is a monomial and c the associate coefficient.\n\n\n\n\n\nseries(f::Function, L::Vector{M}) -> Series{C,M}\n\nCompute the generating series sum_x^α in L f(α) z^α for a function  f mathbbN^n rightarrow C and a sequence L of monomials.\n\n\n\n\n\nseries(w:: Vector{C}, P::Matrix{C}, L::Vector{M}) -> Series{C,M}\n\nCompute the series of the moment sequence _i ω_i P_i^α for α in L.\n\n\n\n\n\nseries(w:: Vector{C}, P::AbstractMatrix, X, d::Int64) -> Series{C,M}\n\nCompute the series of the moment sequence _i ω_i P_i^α for α leq d.\n\n\n\n\n\nseries(p::Polynomial, zeta, X, d::Int64) -> Series\n\nCompute the series of moments p(ζ^α) for α leq d.\n\n\n\n\n\nseries(H::Matrix{C}, L1::Vector{M}, L2::Vector{M}) -> Series{C,M}\n\nCompute the series associated to the Hankel matrix H, with rows (resp. columns) indexed by the array of monomials L1 (resp. L2).\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.scale",
    "page": "Series",
    "title": "MultivariateSeries.scale",
    "category": "function",
    "text": "Scale the moments σ_α by λ^deg(α).\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.scale!",
    "page": "Series",
    "title": "MultivariateSeries.scale!",
    "category": "function",
    "text": "Scale the moments σ_α by λ^deg(α), overwriting σ\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#Base.:*",
    "page": "Series",
    "title": "Base.:*",
    "category": "function",
    "text": "Multiply the elements of L by the variable v\n\n\n\n\n\n *(v::Variable,   σ::Series{C,M}) -> Series{C,M}\n *(m::Monomial,   σ::Series{C,M}) -> Series{C,M}\n *(t::Term,       σ::Series{C,M}) -> Series{C,M}\n *(p::Polynomial, σ::Series{C,M}) -> Series{C,M}\n\nThe dual product (or co-product) where variables are inverted in the polynomial and the monomials with positive exponents are kept in the series.\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariatePolynomials.maxdegree",
    "page": "Series",
    "title": "MultivariatePolynomials.maxdegree",
    "category": "function",
    "text": "maxdegree(σ::Series) -> Int64\n\nMaximal degree of the moments defined in the series σ.\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.hankel",
    "page": "Series",
    "title": "MultivariateSeries.hankel",
    "category": "function",
    "text": "hankel(σ::Series{C,M}, L1::Vector{M}, L2::Vector{M}) -> Array{C,2}\n\nHankel matrix of σ with the rows indexed by the list of polynomials L1 and the columns by L2. The entries are the dot product for σ of the corresponding elements in L1 and L2.\n\nExample\n\njulia> L =[1, x1, x2, x1^2, x1*x2, x2^2]\n\njulia> H = hankel(s,L,L)\n6x6 Array{Float64,2}:\n  4.0   5.0   7.0    5.0  11.0  13.0\n  5.0   5.0  11.0   -1.0  17.0  23.0\n  7.0  11.0  13.0   17.0  23.0  25.0\n  5.0  -1.0  17.0  -31.0  23.0  41.0\n 11.0  17.0  23.0   23.0  41.0  47.0\n 13.0  23.0  25.0   41.0  47.0  49.0\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#MultivariateSeries.hankelbasis",
    "page": "Series",
    "title": "MultivariateSeries.hankelbasis",
    "category": "function",
    "text": "Generate the table of Hankel matrices H_α associated to the monomials x^α generating the space of Hankel matrices indexed by L1 (for the rows) and L2 (for the columns).\n\n\n\n\n\n"
},

{
    "location": "code/1.series/#Series-1",
    "page": "Series",
    "title": "Series",
    "category": "section",
    "text": "Pages = [\"1.series.md\"]MultivariateSeries.SeriesdualLinearAlgebra.dot momentseries scalescale!*maxdegreehankel hankelbasis"
},

{
    "location": "code/2.polynomials/#",
    "page": "Polynomials",
    "title": "Polynomials",
    "category": "page",
    "text": ""
},

{
    "location": "code/2.polynomials/#MultivariateSeries.@ring",
    "page": "Polynomials",
    "title": "MultivariateSeries.@ring",
    "category": "macro",
    "text": "@ring args...\n\nDefines the arguments as variables and output their array.\n\nExample\n\nX = @ring x1 x2\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#MultivariateSeries.monoms",
    "page": "Polynomials",
    "title": "MultivariateSeries.monoms",
    "category": "function",
    "text": "monoms(V, d::Int64) -> Vector{Monomial}\nmonoms(V, rg::UnitRangeInt64) -> Vector{Monomial}\n\nList of all monomials in the variables V up to degree d of from degree d1 to d2, ordered by increasing degree.\n\n\n\n\n\nmonoms(V, d::Int64) -> Vector{Monomial}\n\nList of all monomials in the variables V up to degree d of from degree d1 to d2, ordered by increasing degree.\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#Base.inv",
    "page": "Polynomials",
    "title": "Base.inv",
    "category": "function",
    "text": " inv(m :: Monomial{true})\n\nreturn the inverse monomial with opposite exponents.\n\n\n\n\n\n inv(m :: Monomial)\n\nreturn the inverse monomial with opposite exponents.\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#Base.Math.exponent",
    "page": "Polynomials",
    "title": "Base.Math.exponent",
    "category": "function",
    "text": "exponent(m::Monomial) -> Array{Int64,1}\n\nGet the exponent of a monomial as an array of Int64\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#MultivariateSeries.deg",
    "page": "Polynomials",
    "title": "MultivariateSeries.deg",
    "category": "function",
    "text": "deg(p:Polynomial) -> Int64\n\nDegree of a polynomial\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#MultivariateSeries.sparse_pol",
    "page": "Polynomials",
    "title": "MultivariateSeries.sparse_pol",
    "category": "function",
    "text": "sparse_pol(w, E, X) -> Polynomial{true,C}\n\nCompute the polynomial  ωᵢ X^Ei with coefficients ωᵢ and monomial exponents Ei.\n\n\n\n\n\n"
},

{
    "location": "code/2.polynomials/#Polynomials-1",
    "page": "Polynomials",
    "title": "Polynomials",
    "category": "section",
    "text": "Pages = [\"2.polynomials.md\"]@ring monomsinvexponentdegsparse_pol"
},

{
    "location": "code/3.decompose/#",
    "page": "Decomposition",
    "title": "Decomposition",
    "category": "page",
    "text": ""
},

{
    "location": "code/3.decompose/#MultivariateSeries.decompose",
    "page": "Decomposition",
    "title": "MultivariateSeries.decompose",
    "category": "function",
    "text": "decompose(σ :: Series{C,M}, rkf :: Function)\n\nDecompose the series σ as a weighted sum of exponentials. Return ω, Ξ where\n\nω is the vector of weights,\nΞ is the matrix of frequency points, stored per row.\n\nThe list of monomials of degree leq d-1 over 2 are used to construct the Hankel matrix, where d is the maximal degree of the moments in σ.\n\nThe optional argument rkf is the rank function used to determine the numerical rank from the vector S of singular values. Its default value eps_rkf(1.e-6) determines the rank as the first i s.t. S[i+1]/S[i]< 1.e-6 where S is the vector of singular values.\n\nIf the rank function cst_rkf(r) is used, the SVD is truncated at rank r.\n\n\n\n\n\n"
},

{
    "location": "code/3.decompose/#Decomposition-1",
    "page": "Decomposition",
    "title": "Decomposition",
    "category": "section",
    "text": "Pages = [\"3.decompose.md\"]decompose "
},

]}
