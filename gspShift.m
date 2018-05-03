function [ Y ] = gspShiftAep( signal, N, type, U )
%shift graph signal N times, according to shift operator specified in
%'type'
%returns shifted graph signal
%U is set of eigenvectors, 





    %working on energy perserving shifting
    k = 1:N;
    lambda_e = exp(1i*((-2*pi*(k-1))/N));

    %trying some different things for shift operators
    S = U*diag(lambda_e)*inv(U);
    
    Y = S^N*signal;
    
    
    
end


end
