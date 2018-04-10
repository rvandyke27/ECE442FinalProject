%Extracting edges from gml file graph
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

file = importdata("kangaroo.dat");
A = file.data;
 size(A)
d = zeros((length(A(1,:))),1);
%build degree matrix
for i=1:length(A(:,1))
 for j=1:length(A(1,:))
     if(A(i,j) > 0)
         d(i) = d(i) + 1;
     end

 end
end

D = diag(transpose(d));
size(d)
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

for i = 1:N
    
    F(i, 1) = dot(signal, U(:,i));
end

figure(1)
stem(lambda_vector, F);
title("Graph Fourier Transform for Signal as Degree Matrix");
xlabel('$$\lambda_l$$','Interpreter','Latex');
ylabel('$$\hat{f}(\lambda_l)$$','Interpreter','Latex');

lambda_vector

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

figure(2)
nodes = 1:N;
stem(nodes, f);
title("Recovered Graph Signal from Inverse Fourier Transform");
xlabel("Graph Vertices, Vi");
ylabel("f(Vi)");






