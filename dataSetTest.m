%Testing with additional data sets

[A, D, sig] = getGraphFromGML("karate.gml");
N = length(A);%reconstructed graph signal

L = D - A;

%get eigenvectors (U) and eigenvalues(lambda)
[U, lambda] = eig(L);

F = GFT(U, lambda, sig, N);

f = iGFT(F, U, N);

[x, xs] = SelectionSampling(35, 35, U, N, sig);

[xr, y1] = AggregationSampling(35, U, lambda, N, sig);
