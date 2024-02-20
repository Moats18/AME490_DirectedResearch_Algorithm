function yMin = minY(x, y, Fj, Tj, J, R, A, e)
% 
% uses the method of lagrange multipliers to convert the quadratic 
% energy equation into a linear equation that can be solved 
% relatively easily
%
% inputs:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
% e: e is a vector of the rigid lengths between different y coordinates
%
% Tj: a cell array of the set of all x's within each panel(the jth panel
% corresponds to the jth row)
%
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a 2D array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)

% ***** update y structure***
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% J: the set of all panels
% R: a cell array of all of the rotation matrices for each panel
%
% outputs:
% yMin: y coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   

% Define the matrices that allow the singular vector yi to be converted to
% a list of vectors called y

n = length(y);
lenJ = length(J);

for i = 1:n
    x = zeros(3, 3*n);
    x(1:3, 3*i-2:3*i) = eye(3);
    xM{i} = x;
end 

% Defining the matrix that allows the vector y to be factored out
% Aij{1, 1} is a 3 by 3*n matrix

for j = 1:lenJ
    l = length(Fj{j});
    sum = (1/l)*calcMatrixSum(xM, Fj(j, :));
    for i = 1:n
        Aij{i, j} = xM{i} - sum;
    end
end

% pre-allocating the size of the kMatrix (3n by 3n)
kMatrix = zeros(3*n, 3*n);

for j = 1:lenJ
    for i = length(Fj(j, :))
    f = Fj(j, :);
    k = f(i);
    kMatrix = kMatrix + 2*Aij{k, j}'*Aij{k, j};
    end
end

% pos vectors with respect to the center of the panel for all panels
rij = zeros(3*n, 1, lenJ);
for j = 1:lenJ
[~, rij(:,1,j)] = centerOfPanel(Tj(j, :), x);
end 

% calculation of the b vector
bVector = zeros(3*n, 1);
for j = 1:lenJ
    for i = length(Fj(j, :))
        f = Fj(j, :);
        k = f(i);
        bVector = bVector + rij(k, 1, j)'*R{j}'*Aij{k, j} ; 
    end
end

% B is a 1 by 3*n vector (using the dot product function in MATLAB
% eliminates the need to utilize the tranpose
B = -2*bVector; 

% calculation of the c scalar
cScalar = 0; 
for j = 1:lenJ
    for i = length(Fj(j, :))
        f = Fj(j, :);
        k = f(i);
        cScalar = cScalar + norm(rij(k, 1, j))^2; 
    end
end

% Utilize the rigidity constraint matrix
% write the final solution using the derived form

sol = [kMatrix, A', A, 0] * [B; e];

% take the first 3*n by 1 elements which coorespond to the minimized y
% values

yMin = sol(1:3*n, 1);

end

