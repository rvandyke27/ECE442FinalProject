function [ xr, Vk ] = SelectionSampling( P, K, U, N, signal )
%SelectionSampling with number of observations P, number of eigenvectors
%used K, eigenvectors U, and number of nodes N
%   Returns samples Vk and displays reconstructed signal

C = eye(N);
C = C(1:P,:);
xs = C*signal;

Vk = U(:,1:K);
xr = Vk*transpose(C*Vk)*xs;
figure
stem(xr);
title("Selection Sampling Recovered Graph Signal");

end

