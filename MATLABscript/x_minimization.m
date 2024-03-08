function Xmin = x_minimization(x, y, Fj, Sj, J, R,U,h)
%
%Inputs:
%
% Rigidity Constraints:
% U: U is a matrix that contains the rigidity constraints such that Ux = h
% h: h is a vector of the rigid lenghts between different x coordinates\
%
% Indexing Inputs:
% x: x coordinate 2-D array (3*m by 1 where m is the number of indices)
% y: y coordinate 2-D array (3*m by 1 where m is the number of indices)
% Fj: a cell array of the set of all y's wihtin each panel
% Sj: a cell array of the set of all y's within each panel
% J: the set of all panels
% R: a cell array of all of the rotation matrices for each panel
% 
% Outputs:
%Xmin: x coordinate(2-D array) that minimizes the energy based on given
%constraints

% Define the matrices that allow the singular vector xi to be converted to
% a list of vectors called x

m = length(x);
lenJ = length(J);
mm = length(Fj);

for i = 1:(m/3)
    x = zeros(3,m);
    x(1:3,3*i-2:3*i) = eye(3);
    zM{i} = x;
end

% Defining the matrix that allows the vector x to be factored out
% Gij{1, 1} is a 3 by n matrix

for j = 1:lenJ
    l = length(Fj(:, :, j));
    part_sum = (R{j}/l)*calcMatrixSum(zM, Fj(:, :, j));
    for i = 1:m/3
        Gij{i,j} = R{j}*zM{i} - part_sum;
    end
end

% pos vectors with respect to the center of the panel
dij = zeros(3*mm,1,lenJ);
for j = 1:lenJ
[~, dij(:,:,j)] = centerOfPanel(Sj(:,:,j), y);
end

R{1}
Gij{1, 1}
dij(1:3, 1, 1)
dij(1:3, 1, 1)'*R{1}'*Gij{1, 1}
%calculation of m vector
mvector = zeros(m,1);
for len = 1:lenJ
    for i = length(Fj(:,:,j))
        p = Fj(:,:,j);
        q = p(i);
        mvector = mvector + Gij{q,j}'*dij(3*i-2:3*i, 1, j);
    end
end

%calculation of gmatrix (3m by 3m)
gmatrix = zeros(m,m);
for j = 1:lenJ
    for i=length(Fj(:,:,j))
        p = Fj(:,:,j);
        q = p(i);
        gmatrix = gmatrix + 2*Gij{q,j}'*Gij{q,j};
    end
end

% M is a 1 by 3*m vector (using the dot product function in MATLAB
% eliminates the need to utilize the tranpose
M = -2*mvector; 

%calculating d scalar 
dscalar = 0;
for j = 1:lenJ
    for i=length(Fj(:,:,j))
        p = Fj(:,:,j);
        q = p(i);
        dscalar = dscalar + norm(dij(q,1,j))^2;
    end
end
% Utilize the rigidity constraint matrix
% write the final solution using the derived form
sol = zeros(3*m + length(h),1);

%stores rows and columns of the matrix U
siz = size(U);
%stores the rows of the matrix 
row = siz(1);
big_matrix = zeros(m+row,m+row);
big_matrix(1:m,1:m) = gmatrix;
big_matrix(m+1:m+row,1:m) = U;
big_matrix(1:m,m+1:m+row) = U';
size(M)
size(h)
sol = big_matrix\[M;h];

% taking the first 3*m by 1 element that corresponds to the final minimized
% X value
Xmin = sol(1:m,1);
