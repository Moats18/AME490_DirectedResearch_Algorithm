function yMin = minY(x, y, Fj, Tj, J, R)
% 
% uses the method of lagrange multipliers to convert the quadratic 
% energy equation into a linear equation that can be solved 
% relatively easily
%
% inputs:
%
% need the rigidity contraint
%
% Tj: a cell array of the set of all x's within each panel 

% ***** update x structure ***
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a cell array of the set of all y's within each panel 

% ***** update y structure***
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% J: the set of all opanels
% R: a cell array of all of the rotation matrices for each panel
%
% outputs:
% yMin: 


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

for j = 1:lenJ
    l = length(Fj{j});
    sum = (1/l)*calcMatrixSum(xM, Fj{j});
    for i = 1:n
        Aij{i, j} = xM{i} - sum;
    end
end

% pre-allocating the size of the kMatrix (3n by 3n)
kMatrix = zeros(3*n, 3*n);

for j = 1:lenJ
    for i = length(Fj{j})
    f = Fj{j};
    k = f(i);
    kMatrix = kMatrix + 2*Aij{k, j}'*Aij{k, j};
    end
end

% pos vectors with respect to the center of the panel
for j = 1:lenJ
[~, rij{:, j}] = centerOfPanel(Tj{j}, x);
end 

% calculation of the b vector
bVector = zeros(3, 1, 15);
for j = 1:lenJ
    for i = length(Fj{j})
        f = Fj{j};
        k = f(i);
        bVector = bVector + rij{i, j}'*R{j}'*Aij{k, j} ; 
    end
end

end

