function [ Y ] = gspShiftAc(signal, shifts)
%shift graph signal 'shifts' times by circulant matrix
%returns shifted version of signal

x = transpose(signal);
N = length(x);
Ac = zeros(length(signal));
for i=1:length(Ac(:,1))
 for j=1:length(Ac(1,:))
     if(mod(i - j, N) == 1)
         Ac(i,j) = 1;
     end

 end
end

S = Ac;

Y = S^shifts*x;
    

end