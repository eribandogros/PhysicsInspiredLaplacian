%sample [-1, 1] x [-1, 1]
l = 0.5 %spacing
bound1 = -5;
bound2 = 5;

n = (bound2-bound1)/l;
nV = (n+1)*(n+1);
nE = 2*n*(n+1);

Xarray = (bound1:l:bound2);
%positions organized bottom to top, left to right
pos = zeros(nV,2);
indexV = 1;
for i=1:n+1 %x
    tmp = Xarray(i);
    for j=1:n+1 %y
        pos(indexV,1) = tmp;
        pos(indexV,2) = Xarray(j);
        indexV=indexV+1;
    end
end

X = diag(pos(:,1));
Y = diag(pos(:,2));

%build Hodge Laplacian operators
%S0 = (1/l^2)*eye(nV);
S1 = eye(nE);
D0=zeros(nE,nV);

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
D0;

%Laplacian
LH = D0'*S1*D0;

%Hamiltonian
k = 5; %stiffness
m = 1; %mass
omega = sqrt(k/m)
HX = (1/(l^2))*(.5/m)*LH + .5*(m*omega^2)*X^2;
HY = .5*(m*omega^2)*Y^2;
H = HX + HY;

%solve for eigenvalues
eigs(H,12,'sm')
