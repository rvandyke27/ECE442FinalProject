%ECE 442 Final Project
%Avi Webberman, Rebecca Van Dyke

%Extracting edges from gml file graph
fileName = 'power.gml';
inputfile = fopen(fileName);

A=[];
G=[];
l=0;
k=1;

while 1
    % Get a line from the input file
    tline = fgetl(inputfile);

    % Quit if end of file
    if ~ischar(tline)
        break
    end
    
    nums = regexp(tline,'\d+','match');
    
    if length(nums)
        if l==1
            l=0;
            G(k,2)=str2num(nums{1});  
            k=k+1;
            continue;
        end
        
        G(k,1)=str2num(nums{1});
        l=1;

    else
        l=0;
        continue;
    end

end

%build adjacency matrix

for i=1:length(G)
    A(G(i,1)+1,G(i,2)+1) = 1;
    A(G(i,2)+1,G(i,1)+1) = 1;
end

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
 
%calculate laplacian matrix
L = D - A;

 