function [ f ] = iGFT( F, U, N )
%iGFT computes and displays inverse graph Fourier transform given
% eigenvectors and eigenvalues, graph signal, and N

%reconstructed graph signal
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

end

