function Xmin = x_minimization(x, y, Fj, Uj, J, R)
%
%inputs:
%
% x: x coordinate 2-D array (3*m by 1 where m is the number of indices)
% Fj: a cell array of the set of all y's wihtin each panel
% y: y coordinate 2-D array (3*m by 1 where m is the number of indices)
% J: the set of all panels
% I: the set of all vertices
% R: a cell array of all of the rotation matrices for each panel
% Uj: a cell array of the set of all y's within each panel
% outputs:
%Xmin:

% Define the matrices that allow the singular vector xi to be converted to
% a list of vectors called x

m = length(x);
lenJ = length(J);

for i = 1:m
    x = zeros(3,3*m);
    x(1:3,3*i-2:3*i) = eye(3);
    zM{i} = x;
end

% Defining the matrix that allows the vector x to be factored out

for j = 1:len(J)
    l = length(Fj{j});
    part_sum = (R{j}/l)*x_calcMatrixSum(zM, Fj{j});
    for i = 1:m
        Gij{i,j} = R{j}*zM{i} - part_sum;
    end
end

% pos vectors with respect to the center of the panel
for j = 1:lenJ
[~, dij{:, j}] = x_centerOfPanel(Uj{j}, y);
end

%calculation of m vector
mvector = zeros(3,1,15);
for len = 1:lenJ
    for i = length(Fj{j})
        p = Fj{j};
        q = p{i};
        mvector = mvector + 2*Gij{q,j}'*dij{i,j};
    end
end

%calculation of gmatrix (3m by 3m)
gmatrix = zeros(3*m,3*m);
for j = 1:lenJ
    for i=length(Fj{j})
        p = Fj{j};
        q = p{i};
        gmatrix = gmatrix + 2*Gij{q,j}*Gij{q,j};
    end
end
end


