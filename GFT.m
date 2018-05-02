function [ F ] = GFT( U, lambda, sig, N )
%GFT computes and displays graph Fourier transform given eigenvectors and eigenvalues,
%graph signal, and N

F = zeros(N, 1);
lambda_vector = zeros(N,1);

for i = 1:N
    lambda_vector(i) = lambda(i,i);
end

for i = 1:N  
    F(i, 1) = dot(transpose(sig), U(:,i));
end

figure(1)
stem(lambda_vector(1:length(lambda_vector)), F(1:length(F)));
title("Graph Fourier Transform for Signal as Degree Matrix");
xlabel('$$\lambda_l$$','Interpreter','Latex');
ylabel('$$\hat{f}(\lambda_l)$$','Interpreter','Latex');

%disp(lambda_vector)

end

