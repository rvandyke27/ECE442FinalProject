function [ xr, y1 ] = AggregationSampling( K, U, lambda, N, signal )
%AggregationSampling with number of eigenvectors used K, eigenvectors U, 
%eigenvalues lambda, number of nodes N, and graph signal
%   Returns samples y1 and displays reconstructed signal

E_kT = [eye(K) zeros(K,N-K)];
C = E_kT;
signal = transpose(signal);

y1 = zeros(K, 1);
for i=1:K
    shifted = gspShiftAep(signal, i, U); 
    y1(i) = shifted(1); 
end

y1_hat = C*y1;

Vk = U(:,1:K);
u1_hat = transpose(Vk)*transpose(E_kT(1,:));
disp(size(u1_hat));

psi = zeros(N);
for l = 1:N
    for k = 1:N
        psi(l,k) = lambda(k)^(l-1);
    end
end

disp(size(C*transpose(psi)*transpose(E_kT)));
disp(size(diag(u1_hat)));
xr = Vk / (diag(u1_hat))*((C*transpose(psi)*transpose(E_kT)))\y1_hat;

figure
stem(xr);
title("Aggregation Sampling Recovered Graph Signal");

end

