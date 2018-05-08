excelfile = 'gavett_long_02_edit2.csv';
timestamp = transpose(csvread(excelfile,1,0,[1 0 1616 0]));
x = transpose(csvread(excelfile,1,1,[1 1 1616 1]));
y = transpose(csvread(excelfile,1,2,[1 2 1616 2]));
rssi = transpose(csvread(excelfile,1,4,[1 4 1616 4]));
ping = transpose(csvread(excelfile,1,5,[1 5 1616 5]));
% Set RSSI values low to -100
for i=1:length(rssi)
    if rssi(i) == 0
        rssi(i) = -100;
    end
end
sample = 1:2:1616;
rssi_samples = rssi(sample);
ping_samples = ping(sample);
nodes = 1:length(rssi_samples);
coord = [x(sample); y(sample)];

A = zeros(length(nodes));

for i=1:length(nodes)
    for j=1:length(nodes)
        
        distance = sqrt((coord(1,i) - coord(1,j))^2 + (coord(2,i) - coord(2,j))^2);
        if(distance > 0)
            A(i,j) = 1/distance;
        end
        

    end
end


%unweighted adjacency matrix
Auw = zeros(size(A));

 %size(A)
d = zeros((length(A(1,:))),1);
deg = zeros(size(A));
%build degree matrix
for i=1:length(A(:,1))
 for j=1:length(A(1,:))
     if(A(i,j) > 0)
         d(i) = d(i) + A(i,j);
         Auw(i,j) = 1;
     end

 end
end

D = diag(transpose(d));
size(d);

%calculate laplacian matrix
L = D - A;

%get eigenvectors (U) and eigenvalues(lambda)
[U, lambda] = eig(L);

size(A);
size(U);
N = length(d);
F = zeros(N, 1);
Fping = zeros(N,1);
lambda_vector = zeros(N,1);


for i = 1:N
    lambda_vector(i) = lambda(i,i);
end

%10 denotes a female kangaroo and -10 denotes a male kangaroo
signal = rssi_samples; 
signal_ping = ping_samples;
%signal = d;

for i = 1:N  
    F(i, 1) = dot(signal, U(:,i));
    Fping(i,1) = dot(signal_ping, U(:,i));
end

figure(1)
stem(lambda_vector(1:length(lambda_vector)), F(1:length(F)));
title("Graph Fourier Transform for Signal as Degree Matrix");
xlabel('$$\lambda_l$$','Interpreter','Latex');
ylabel('$$\hat{f}(\lambda_l)$$','Interpreter','Latex');

disp(lambda_vector);
U(:,1);

%reconstructed graph signal
f = zeros(N,1);
fping = zeros(N,1);

%compute inverse fourier transform
for i = 1:N
    sum = 0;
    sum2 = 0;
    for j = 1:N
        sum = sum + U(i,j)*F(j);
        sum2 = sum2 + U(i,j)*Fping(j);
    end
    f(i) = sum;
    fping(i) = sum2;
end

%plot recovered plot 
figure(2)
nodes = 1:N;
stem(nodes, f);
title("Recovered Graph Signal from Inverse Fourier Transform");
xlabel("Graph Vertices, Vi");
ylabel("f(Vi)");



 
Ac = zeros(size(A));
for i=1:length(Ac(:,1))
 for j=1:length(Ac(1,:))
     if(mod(i - j, N) == 1)
         Ac(i,j) = 1;
     end

 end
end


x = transpose(signal);
lambda_max = max(lambda_vector);

%working on energy perserving shifting
k = 1:N;
lambda_e = exp(1i*((-2*pi*(k-1))/N));

%trying some different things for shift operators
Aphi = U*diag(lambda_e)*inv(U);
S1 = Ac;
S2 = Aphi;
S3 = L;
S4 = A;

Y1 = S1^20*x;
Y2 = (S2)^20*x;
Y3 = S3*x;
Y4 = S4*x;

% Yh = H*x;

%Plotting shifted graph signal
figure(3)
hold on
stem(nodes, Y1, 'r');
title("Shifted");
hold off

%Graph filters

%polynomial filter stuff
filt = N;
cutoff = 10;
h = zeros(filt, 1);
h(length(h)-cutoff:length(h)) = 1/(cutoff+1);
%h(1:cutoff) = 1/(cutoff+1);
hp = transpose(transpose(h));
Hp = zeros(size(S1));

% for l = 1:length(hp)
%     Hp = Hp + (hp(l)*S2^(l-1));
% end
% 
% Yhp = Hp*x;
% figure(4)
% stem(nodes, Yhp);
% title("Filtering (Polynomial)");

%frequency response of filter
vmonde = zeros(length(hp), N);

for l = 1:length(hp)
    for k = 1:N
        vmonde(l,k) = lambda_e(k)^(l-1);
    end
end



hf = transpose(vmonde)*hp;

Hf = U*diag(hf)*inv(U);

Yhf = Hf*x;

figure(5)
stem(nodes, Yhf);
title("Signal after High Pass Filter")


%Sampling


%Selection Sampling
P = 800;
K = 600;
C = eye(N);
C = C(1:P,:);
xs = C*x;

Vk = U(:,1:K);
xr = Vk*transpose(C*Vk)*xs;
figure(6)
stem(nodes, xr);
title("Recovered Graph Signal from Selection Sampling (using 100 Eigenvectors)");

%Aggregation Sampling

figure(7)
stem(nodes, x)
title("Initial Graph Signal (RSSI)");

figure(8)
stem(lambda_vector, Fping);


%-----------------------------------------------------------
%-----------------------------------------------------------
function [ xr, y1 ] = AggregationSampling( K, U, lambda, N, signal )
%AggregationSampling with number of eigenvectors used K, eigenvectors U, 
%eigenvalues lambda, number of nodes N, and graph signal
%   Returns samples y1 and displays reconstructed signal

E_kT = [eye(K) zeros(K,N-K)];
C = E_kT;
signal = transpose(signal);


y1 = zeros(N, N);

for i=1:N
    shifted = gspShiftAep(signal, i, U); 

    y1(i,:) = transpose(shifted);
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

% disp(size(C*transpose(psi)*transpose(E_kT)));
% disp(size(diag(u1_hat)));
% xr = Vk / (diag(u1_hat))*((C*transpose(psi)*transpose(E_kT)))\y1_hat;
xr = Vk*pinv(diag(u1_hat))*pinv(C*transpose(psi)*transpose(E_kT))*y1_hat;
size(xr)
% figure
% stem(xr);
% title("Aggregation Sampling Recovered Graph Signal");

end

