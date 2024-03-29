function [ A, D, sig] = getGraphFromGML( fileName )
%getGraphFromGML Returns A, adjacency matrix of graph from GML file
%d, degree matrix
%and sig, graph signal values computed by in-degree

% Extracting edges from gml file graph
inputfile = fopen(fileName);

l=0;
k=1;
G = [];

while 1
    % Get a line from the input file
    tline = fgetl(inputfile);
    % Quit if end of file
    if ~ischar(tline)
        break
    end
    
    nums = regexp(tline,'\d+','match');

 
    if ~isempty(nums)

        if l==1
            l=0;
            G(k,2)=str2double(nums{1});  
            k=k+1;
            continue;
        end
        G(k,1)=str2double(nums{1});
        l=1;
    else
        l=0;
        continue;
    end
end

%build adjacency matrix
A = [];
for i=1:length(G)
    A(G(i,1)+1,G(i,2)+1) = 1;
    A(G(i,2)+1,G(i,1)+1) = 1;
end

Auw = zeros(size(A));

%size(A)
d = zeros((length(A(1,:))),1);

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

%define graph signal (popularity)
sig = zeros(size(A, 1), 1);
for i=1:length(sig)
   sig(i, 1) = sum(A(i, :)); 
end
 
end

