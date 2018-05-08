function [ Yf ] = gspFilter(signal, type, cutoff, U)
%shift graph signal 'shifts' times by circulant matrix
%returns shifted version of signal
filt = length(signal);
h = zeros(filt, 1);
N = length(signal);
x = transpose(signal);
if(type == 'lpf')

    h(1:cutoff) = 1/(cutoff+1);
elseif(type == 'hpf')
    h(length(h)-cutoff:length(h)) = 1/(cutoff+1);
else
    error("Please enter valid filter type: 'lpf' or 'hpf'");
end

%working on energy perserving shifting
k = 1:N;
lambda_e = exp(1i*((-2*pi*(k-1))/N));

vmonde = zeros(length(h), N);

for l = 1:length(h)
    for k = 1:N
        vmonde(l,k) = lambda_e(k)^(l-1);
    end
end



hf = transpose(vmonde)*h;

Hf = U*diag(hf)*inv(U);

Yf = Hf*x;

    

end