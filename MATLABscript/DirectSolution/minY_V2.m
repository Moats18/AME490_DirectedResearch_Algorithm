function yMin = minY_V2(x, y, Fj, Tj, J, R, A)
% 
% directly solves the linear equation to calculate the perturbation of the
% initial y 
%
% Inputs
%
% Rigidity Constraints:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
%
% Indexing Inputs:
% Tj: a cell array of the set of all x's within each panel(the jth panel
% corresponds to the jth row)
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a 2D array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% J: the set of all panels
% R: a cell array of all of the rotation matrices for each panel
%
% Outputs
% yMin: y coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   
% Define the matrices that allow the 3D coord vector yi to be converted to
% a list of coord vectors called y

% creating the k matrix
n = length(y);
lenJ = length(J);
nn = length(Fj);

for i = 1:(n/3)
    x1 = zeros(3, n);
    x1(1:3, 3*i-2:3*i) = eye(3);
    xM{i} = x1;
end 

% Defining the matrix that allows the vector y to be factored out
% Aij{1, 1} is a 3 by n matrix

for j = 1:lenJ
    l = length(Fj(:, :, j));
    sum = (1/l)*calcMatrixSum(xM, Fj(:, :, j));
    for i = 1:length(Fj(:, :, j))
        k = Fj(:, i, j);
        Aij{i, j} = xM{k} - sum;
    end
end

% pre-allocating the size of the kMatrix (3n by 3n)
kMatrix = zeros(n, n);

for j = 1:lenJ
    for i = 1:length(Fj(:, :, j))
    kMatrix = kMatrix + 2*Aij{i, j}'*Aij{i, j};
    end
end

%determining the N matrix
nMatrix = null(A);

% pos vectors with respect to the center of the panel for all panels
rij = zeros(3*nn, 1, lenJ);

for j = 1:lenJ
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

% calculation of the negative b vector
bVector = zeros(n, 1);
for j = 1:lenJ
    for i = 1:length(Fj(:, :, j))
        %(rij(3*i-2:3*i, 1, j)'*R{j}'*Aij{i, j})'
        bVector = bVector + (rij(3*i-2:3*i, 1, j)'*R{j}'*Aij{i, j})'; 
    end
end

% B is a 1 by 3*n vector 
B = -2*bVector; 

% calculation of bTilde
bTilde = -1*(nMatrix'*kMatrix*y + nMatrix'*B); % removed mutlipe of 2 (for some reason)

% determining the perturbation method
bTilde(abs(bTilde)<1e-3)=0;

yTilde = pinv(nMatrix'*kMatrix*nMatrix)*bTilde;

yMin = y + nMatrix*yTilde;

end

