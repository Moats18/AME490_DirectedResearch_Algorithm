function xMin = minX_V2(x, y, Fj, Tj, J, R, U)
%
% directly solves the linear equation to calculate the perturbation of the
% initial x
%
% Inputs:
%
% Rigidity Constraints:
% U: U is a matrix that contains the rigidity constraints such that Ux = h
%
% Indexing Inputs:
% x: x coordinate 2-D array (2*n by 1 where n is the number of vertices)
% y: y coordinate 2-D array (3*n by 1 where n is the number of vertices)
% Tj: a cell array of the set of all x's within each panel (the jth panel
% corresponds to the jth row)
% Fj: a cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% J: the labeling set of all panels
% R: a cell array of all of the rotation matrices for each panel
% 
% Outputs:
% Xmin: x coordinate(2-D array) that minimizes the energy based on given
% constraints

% Initialize relevant parameters
n = length(x); % number of vertices * 2
lenJ = length(J); % number of panels
chi = cell(n/2,1); % cell array mapping column X to each vertex x

for k = 1:(n/2)
    chi_k = zeros(3, n);
    chi_k(1:3, 2*k-1:2*k) = [1  0;
                             0  1;
                             0  0];
    chi{k} = chi_k;
end

% Defining the matrix that allows the vector x to be factored out
% Sij{i, j} is a 3 by n matrix, has i rows where i is the maximum number of
% vertices in a panel for the given configuration, and j is the number of
% panels

RjSij = cell(max(cellfun('size', Tj, 2)), lenJ);
for j = 1:lenJ
    l = length(Tj{j});
    sum = R{j}*(1/l)*calcMatrixSum_x(chi, Tj{j}); 
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k
        RjSij{i,j} = (R{j}*chi{k}) - sum; 
    end
end

% pre-allocating the size of A (2n by 2n)
A = zeros(n,n);

for j = 1:lenJ
    for i= 1:length(Fj{j})
        A = A + 2*RjSij{i,j}'*RjSij{i,j};
    end
end

% determining the N matrix
N = null(U);

% pos vectors with respect to the center of the panel for all panels
rij = cell(lenJ); % j cells with the j-th cell containing 
for j = 1:lenJ
    [~, rij{j}] = centerOfPanel3D(Fj{j}, y);
end

% calculation of a
a = zeros(n,1);
for j = 1:lenJ
    for i = 1:length(Fj{j})
       a = a + (RjSij{i,j}'*rij{j}(3*i-2:3*i));
    end
end
a = 2*a; 

% calculation of aTilde
aTilde = N'*a - N'*A*x ; 

% determining the perturbation method
aTilde(abs(aTilde)<1e-5)=0;

xTilde = pinv(N'*A*N)*aTilde;

xMin = x + N*xTilde;

end