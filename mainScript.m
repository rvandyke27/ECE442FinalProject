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

lambda_vector
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


% %random test filter
%  h = transpose([1 1 0.5]);
% %random test filter
%  H = h(1)*S^0 + h(2)*S^1 + h(3)*S^2;

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
stem(nodes, Y2, 'r');
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

vmonde = zeros(length(h), N);

for l = 1:length(h)
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






