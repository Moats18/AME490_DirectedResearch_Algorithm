 %
% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Antoine Moats, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 10/20/24
%
% The initial configuration consists of two panels that share common
% edges. All of the nine vertices of the panels can be determined based on two user 
% specified points. Using the symmetry constraints imposed on six points 
% results in nine fully defined points

% Initial x-values
x1 = [0; 0];
x2 = [1; 0];
x3 = [2; 0];
x4 = [2; 1];
x5 = [1; 1];
x6 = [0; 1];

% x rigidity constraint matrix
U = [-eye(2,2), zeros(2,2), eye(2,2), zeros(2,2), zeros(2,2), zeros(2,2);
     zeros(2,2), zeros(2,2), zeros(2,2), eye(2,2), zeros(2,2), -eye(2,2);
     -eye(2,2), zeros(2,2), eye(2,2), zeros(2,2), zeros(2,2), zeros(2,2);
     -eye(2,2), zeros(2,2), zeros(2,2), zeros(2,2), zeros(2,2), eye(2,2);
     zeros(2,2), zeros(2,2), -eye(2,2), eye(2,2), zeros(2,2), zeros(2,2);
     zeros(2,2), -eye(2,2), zeros(2,2), zeros(2,2), eye(2,2), zeros(2,2)];

% setting the x and y constraint matrix to be the same i.e. stating that
% they have the same overall shape
A = [-eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
      zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3), zeros(3,3), -eye(3,3);
     -eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     -eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3);
      zeros(3,3), zeros(3,3), -eye(3,3), eye(3,3), zeros(3,3), zeros(3,3);
      zeros(3,3), -eye(3,3), zeros(3,3), zeros(3,3), eye(3,3), zeros(3,3)];


% final x vector
x = [x1; x2; x3; x4; x5; x6];

% Initial y-values
y1 = [0; 0; 0];
y2 = [0.75; 0; 1];
y3 = [1.5; 0; 0];
y4 = [1.5; 1; 0];
y5 = [0.75; 1; 1];
y6 = [0; 1; 0];

% final y vector
y = [y1; y2; y3; y4; y5; y6];

% populate the vector numbering all of the panels
J = [1, 2];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
% a cell array of the set of all x's within each panel 
T1 = F1;
T2 = F2;

% a 3-d array that contains the index set of x coordinates within a given
% panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;

Fj = Tj;

% Initial R 
for j = 1:length(J)
R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for minimization
tol = 10^(-15);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Folding");