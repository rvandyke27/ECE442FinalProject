%Testing with additional data sets

[A, D, sig] = getGraphFromGML("karate.gml");

N = length(A);%reconstructed graph signal
nodes = 1:N;
L = D - A;

%get eigenvectors (U) and eigenvalues(lambda)
[U, lambda] = eig(L);

% F = GFT(U, lambda, sig, N);
% 
% f = iGFT(F, U, N);
% 
%  [x, xs] = SelectionSampling(30, 30, U, N, sig);

[xrec, y1] = AggregationSampling(2, U, lambda, N, sig);
size(xrec)
figure()
stem(nodes, sig);
figure()
stem(nodes, real(xrec(1,:)));