function xMin = minX_V2(x, y, Fj, Sj, J, R, U)
%
% directly solves the linear equation to calculate the perturbation of the
% initial x
%
% Inputs:
%
% Rigidity Constraints:
% U: U is a matrix that contains the rigidity constraints such that Ux = h
% h: h is a vector of the rigid lenghts between different x coordinates\
%
% Indexing Inputs:
% x: x coordinate 2-D array (2*m by 1 where m is the number of indices)
% y: y coordinate 2-D array (3*m by 1 where m is the number of indices)
% Fj: a cell array of the set of all y's wihtin each panel
% Sj: a cell array of the set of all y's within each panel
% J: the set of all panels
% R: a cell array of all of the rotation matrices for each panel
% 
% Outputs:
% Xmin: x coordinate(2-D array) that minimizes the energy based on given
% constraints

%creating the g matrix 
m = length(x);
lenJ = length(J);
mm = length(Fj);

for i = 1:(m/2)
    x1= zeros(3, m);
    x1(1:3, 2*i-1:2*i) = [1  0;
                         0  1;
                         0  0];
    zM{i} = x1;
end

% Defining the matrix that allows the vector x to be factored out
for j = 1:lenJ
    l = length(Fj(:, :, j));
    part_sum = R{j}*(1/l)*calcMatrixSum_x(zM, Fj(:, :, j)); 
    for i = 1:length(Fj(:,:,j))
        k=Fj(:,i,j);
        Gij{i,j} = (R{j}*zM{k}) - part_sum; 
    end
end

%calculation of gmatrix
gmatrix = zeros(m,m);

for j = 1:lenJ
    for i= 1:length(Fj(:,:,j))
        gmatrix = gmatrix + 2*Gij{i,j}'*Gij{i,j};
    end
end

%determining the N matrix
nMatrix = null(U);

% pos vectors with respect to the center of the panel for all panels
dij = zeros(3*mm, 1, lenJ);

for j = 1:lenJ
[~, dij(:,:,j)] = centerOfPanel3D(Sj(:,:,j), y);
end


%calculation of m vector
mvector = zeros(m,1);
for j = 1:lenJ
    for i = 1:length(Fj(:,:,j))
       mvector = mvector +  (Gij{i,j}'*dij(3*i-2:3*i, 1, j));
    end
end

M = -2*mvector; 

% calculation of mTilde
mTilde = -1 * (nMatrix'*gmatrix*x + nMatrix'*M); 

% determining the perturbation method
mTilde(abs(mTilde)<1e-5)=0;

xTilde = pinv(nMatrix'*gmatrix*nMatrix)*mTilde;

xMin = x + nMatrix*xTilde;

end
