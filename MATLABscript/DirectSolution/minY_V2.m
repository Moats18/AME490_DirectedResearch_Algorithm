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
% x: x coordinate 2-D array (2*n by 1 where n is the number of vertices)
% y: y coordinate 2-D array (3*n by 1 where n is the number of vertices)
% Tj: a cell array of the set of all x's within each panel (the jth panel
% corresponds to the jth row)
% Fj: a column cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% J: the labeling set of all panels
% R: a cell array of all of the rotation matrices for each panel
%
% Outputs
% yMin: y coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   

% Initialize relevant parameters
n = length(y); % number of vertices * 3
lenJ = length(J); % number of panels
chi = cell(n/3,1); % cell array mapping column Y to each vertex y

for k = 1:(n/3)
    chi_k = zeros(3, n);
    chi_k(1:3, 3*k-2:3*k) = eye(3); % mapping is identity for the 3 desired coordinate values
    chi{k} = chi_k;
end 

% Defining the matrix that allows the vector y to be factored out
% Sij{i, j} is a 3 by n matrix, has i rows where i is the maximum number of
% vertices in a panel for the given configuration, and j is the number of
% panels

Sij = cell(max(cellfun('size', Fj, 2)), lenJ);
for j = 1:lenJ
    l = length(Fj{j}); % number of vertices in the j-th panel
    sum = (1/l)*calcMatrixSum_y(chi, Fj{j});
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k
        Sij{i, j} = chi{k} - sum;
    end
end

% pre-allocating the size of B (3n by 3n)
B = zeros(n, n);

for j = 1:lenJ
    for i = 1:length(Fj{j})
        B = B + 2*Sij{i, j}'*Sij{i, j};
    end
end

% determining the N matrix
N = null(A);

% pos vectors with respect to the center of the panel for all panels
rij = cell(lenJ); % j cells with the j-th cell containing 
for j = 1:lenJ
    [~, rij{j}] = centerOfPanel2D(Tj{j}, x);
end 

% calculation of b
b = zeros(n, 1);
for j = 1:lenJ
    for i = 1:length(Fj{j})
        r_temp = rij{j}(2*i-1:2*i);
        rij1 = r_temp(1);
        rij2 = r_temp(2);
        b = b + ([rij1 rij2 0]*R{j}'*Sij{i, j})'; 
    end
end
b = 2*b;

% calculation of bTilde
bTilde = N'*b - N'*B*y;

% determining the perturbation method
bTilde(abs(bTilde)<1e-5)=0;

yTilde = pinv(N'*B*N)*bTilde;

yMin = y + N*yTilde;

end