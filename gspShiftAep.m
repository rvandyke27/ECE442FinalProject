function [ Y ] = gspShiftAep( signal, N, U )
%shift graph signal N times by energy preserving adjacency matrix
%U represents set of eigenvectors
%returns shifted graph signal


x = transpose(signal);
%working on energy perserving shifting
k = 1:N;
lambda_e = exp(1i*((-2*pi*(k-1))/N));

%trying some different things for shift operators
S = U*diag(lambda_e)*inv(U);
    
Y = S^N*x;
    
    
    
end



