l = 0.2; %spacing
bound1 = -10;
bound2 = 10;
%smallest spacing recommended is 0.05 for [-5,5] 

n = (bound2-bound1)/l;
nV = (n+1)*(n+1);
nE = 2*n*(n+1);

%atom positions (input)
nA = 2;
A = zeros(nA,2);
A(1,:) = [3 2];
A(2,:) = [-2 -3];

%grid positions bottom to top, left to right
Vdiff = sparse(nV,1);
Vsum = 0;
Xarray = (bound1:l:bound2);
pos = zeros(nV,2);
indexV = 1;
for i=1:n+1 %x
    tmp = Xarray(i);
    for j=1:n+1 %y
        pos(indexV,1) = Xarray(i);
        pos(indexV,2) = Xarray(j);
        
        %potential
        d = zeros(nA,1);
        for k=1:nA
            d(k) = pdist([A(k,:);pos(indexV,:)], 'Euclidean');
            decay_rate = 1/2;
            Vsum = Vsum - (d(k))^decay_rate;
            %Vsum = Vsum + d(k)^decay_rate;
        end

        Vdiff(indexV,1) = Vsum+8;
        
        indexV=indexV+1;
        Vsum = 0;
    end
end
X = diag(pos(:,1));
Y = diag(pos(:,2));

%visualize grid and Vdiff
scatter(pos(:,1),pos(:,2),120,Vdiff,'filled');
colorbar
ylabel('y')
xlabel('x')
title('2D potential for 2 atoms')
xticks(bound1:10*l:bound2); yticks(bound1:10*l:bound2); ax = gca; ax.FontSize = 20;
text(3,2,"\leftarrow atom")
text(-2,-3,"\leftarrow atom")
grid on

%build ops
S1 = speye(nE);
D0=sparse(nE,nV);

edgeIndex = 1;
vertIndex = 0;
for i=1:n+1
    %y edges
    for j=1:n
        D0(edgeIndex,vertIndex+j) = -1;
        D0(edgeIndex,vertIndex+j+1) = 1;
        edgeIndex = edgeIndex+1;
    end
    %x edges
    if i<n+1
        for j=1:n+1
            D0(edgeIndex,vertIndex+j) = -1;
            D0(edgeIndex,vertIndex+j+n+1) = 1;
            edgeIndex = edgeIndex+1;
        end
    end
    vertIndex = vertIndex+n+1;
end

L = D0'*S1*D0;
P = diag(Vdiff(:,1));
H = L+P;

eigs(H,speye(nV),10,'smallestreal')