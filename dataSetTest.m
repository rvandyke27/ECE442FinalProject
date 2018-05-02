%Testing with additional data sets

[A, D, sig] = getGraphFromGML("karate.gml");
N = length(A);%reconstructed graph signal
f = zeros(N,1);

%compute inverse fourier transform
for i = 1:N
    sum = 0;
    for j = 1:N
        sum = sum + U(i,j)*F(j);
    end
    f(i) = sum;
end

%plot recovered plot 
figure(2)
nodes = 1:N;
stem(nodes, f);
title("Recovered Graph Signal from Inverse Fourier Transform");
xlabel("Graph Vertices, Vi");
ylabel("f(Vi)");

L = D - A;

%get eigenvectors (U) and eigenvalues(lambda)
[U, lambda] = eig(L);

F = GFT(U, lambda, sig, N);

f = iGFT(F, U, N);
