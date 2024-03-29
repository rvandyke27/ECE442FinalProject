% % % Extracting edges from gml file graph
% fileName = 'football.gml';
% inputfile = fopen(fileName);
% A=[];
% G=[];
% l=0;
% k=1;
% while 1
%       % Get a line from the input file
%       tline = fgetl(inputfile);
%       % Quit if end of file
%       if ~ischar(tline)
%           break
%       end
%       nums = regexp(tline,'\d+','match');
%       if length(nums)
%           if l==1
%               l=0;
%               G(k,2)=str2num(nums{1});  
%               k=k+1;
%               continue;
%           end
%           G(k,1)=str2num(nums{1});
%           l=1;
%       else
%           l=0;
%           continue;
%       end
% end
% 
% 
%  %build adjacency matrix
%  
%  for i=1:length(G)
%      A(G(i,1)+1,G(i,2)+1) = 1;
%      A(G(i,2)+1,G(i,1)+1) = 1;
%  end
%  

%use kangaroo data set
file = importdata("kangaroo.dat");
A = file.data;

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
lambda_vector = zeros(N,1);


for i = 1:N
    lambda_vector(i) = lambda(i,i);
end

%10 denotes a female kangaroo and -10 denotes a male kangaroo
signal = [10 10 10 10 10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10];
%signal = d;

for i = 1:N  
    F(i, 1) = dot(signal, U(:,i));
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



 
Ac = zeros(size(A));
for i=1:length(Ac(:,1))
 for j=1:length(Ac(1,:))
     if(mod(i - j, N) == 1)
         Ac(i,j) = 1;
     end

 end
end


%random test filter
 h = transpose([1 1 0.5]);
%random test filter
%H = h(1)*S^0 + h(2)*S^1 + h(3)*S^2;

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

Y1 = S1*x;
Y2 = (S2)^N*x;
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
hp = transpose([1 1 1]);
Hp = zeros(size(S1));

for l = 1:length(hp)
    Hp = Hp + (hp(l)*S2^(l-1));
end

Yhp = Hp*x;
figure(4)
stem(nodes, Yhp);
title("Filtering (Polynomial)");

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
title("Filtering (Frequency)")


%Sampling


%Selection Sampling
P = 16;
K = 15;
C = eye(N);
C = C(1:P,:);
xs = C*x;

Vk = U(:,1:K);
xr = Vk*transpose(C*Vk)*xs;
figure(6)
stem(xr);
title("Recovered Graph Signal");

%Aggregation Sampling


%-----------------------------------------------------------
%-----------------------------------------------------------

excelfile = 'gavett_long_02_edit2.csv';
timestamp = transpose(csvread(excelfile,1,0,[1 0 1616 0]));
x = transpose(csvread(excelfile,1,1,[1 1 1616 1]));
y = transpose(csvread(excelfile,1,2,[1 2 1616 2]));
rssi = transpose(csvread(excelfile,1,4,[1 4 1616 4]));
% Set RSSI values low to -100
for i=1:length(rssi)
    if rssi(i) == 0
        rssi(i) = -100;
    end
end
sample = 100:2:1616;
rssi_samples = rssi(sample);
Nodes = 1:length(rssi_samples);
coord = [x(sample); y(sample)];
Ncmr = length(Nodes);
A_cmr = zeros(length(Nodes));
d_cmr = zeros(length(Nodes),1);
for i=1:length(Nodes)
    for j=1:length(Nodes)
        
        distance = sqrt((coord(1,i) - coord(1,j))^2 + (coord(2,i) - coord(2,j))^2);
        if(distance > 0)
            A_cmr(i,j) = 1/sqrt(distance);
        end
        

    end
end

for i=1:length(A_cmr(:,1))
 for j=1:length(A_cmr(1,:))
     if(A_cmr(i,j) > 0)
         d_cmr(i) = d_cmr(i) + A_cmr(i,j);    
     end

 end
end

Lcmr = diag(d_cmr) - A_cmr;

%get eigenvectors (U) and eigenvalues(lambda)
[Ucmr, lambdacmr] = eig(Lcmr);

Fcmr = zeros(Ncmr, 1);
lambda_vectorcmr = zeros(Ncmr,1);


for i = 1:N
    lambda_vectorcmr(i) = lambdacmr(i,i);
end


for i = 1:N  
    Fcmr(i, 1) = dot(rssi_samples, Ucmr(:,i));
end

figure(7)
stem(lambda_vectorcmr(1:length(lambda_vectorcmr)), Fcmr(1:length(Fcmr)));
title("Graph Fourier Transform for RSSI Dataset");
xlabel('$$\lambda_l$$','Interpreter','Latex');
ylabel('$$\hat{f}(\lambda_l)$$','Interpreter','Latex');

%reconstructed graph signal
fcmr = zeros(Ncmr,1);

%compute inverse fourier transform
for i = 1:Ncmr
    sum = 0;
    for j = 1:Ncmr
        sum = sum + Ucmr(i,j)*Fcmr(j);
    end
    fcmr(i) = sum;
end

%plot recovered plot 
figure(8)
stem(Nodes, fcmr(:,1));
title("Recovered Graph Signal from Inverse Fourier Transform");
xlabel("Graph Vertices, Vi");
ylabel("f(Vi)");

figure(9)
stem(Nodes, rssi_samples(1,:))



