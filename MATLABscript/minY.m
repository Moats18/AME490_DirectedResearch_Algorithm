function yMin = minY(x, y, Fj, Tj, J, R, A, e)
% 
% uses the method of lagrange multipliers to convert the quadratic 
% energy equation into a linear equation that can be solved 
% relatively easily
%
% Inputs
%
% Rigidity Constraints:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
% e: e is a vector of the rigid lengths between different y coordinates
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
    for i = 1:n/3
        Aij{i, j} = xM{i} - sum;
    end
end

% pre-allocating the size of the kMatrix (3n by 3n)
kMatrix = zeros(n, n);

for j = 1:lenJ
    for i = length(Fj(:, :, j))
    f = Fj(:, :, j);
    k = f(i);
    kMatrix = kMatrix + 2*Aij{k, j}'*Aij{k, j};
    end
end

% pos vectors with respect to the center of the panel for all panels
rij = zeros(3*nn, 1, lenJ);

for j = 1:lenJ
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

R{1}
Aij{1, 1}
rij(1:3, 1, 1)
rij(1:3, 1, 1)'*R{1}'*Aij{1, 1}
% calculation of the b vector
bVector = zeros(n, 1);
for j = 1:lenJ
    for i = length(Fj(:, :, j))
        f = Fj(:, :, j);
        k = f(i);
        bVector = bVector + (rij(3*i-2:3*i, 1, j)'*R{j}'*Aij{k, j})'; 
    end
end
size(bVector)

% B is a 1 by 3*n vector (using the dot product function in MATLAB
% eliminates the need to utilize the tranpose
B = -2*bVector; 

% calculation of the c scalar
cScalar = 0; 
for j = 1:lenJ
    for i = length(Fj(:, :, j))
        f = Fj(:, :, j);
        k = f(i);
        cScalar = cScalar + norm(rij(k, 1, j))^2; 
    end
end

% Utilize the rigidity constraint matrix
% write the final solution using the derived form
sol = zeros(3*n + length(e), 1);

% stores the rows and columns of matrix A
sz = size(A);
%stores the rows of the matrix 
rows = sz(1);
bigMatrix = zeros(n+rows, n+rows);
bigMatrix(1:n, 1:n) = kMatrix;
bigMatrix(n+1:n+rows, 1:n) = A;
bigMatrix(1:n,n+1:n+rows) = A';
size(B)
size(e)
sol = bigMatrix\[B; e];

% take the first 3*n by 1 elements which coorespond to the minimized y
% values
yMin = sol(1:n, 1);

end

