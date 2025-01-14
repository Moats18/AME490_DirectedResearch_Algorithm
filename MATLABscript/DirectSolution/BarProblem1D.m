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
% The initial configuration consists of two lines that share common
% a common hinge point in the center. In the reference configuration, the hinge is fixed
% resulting in a straight line. In the deformed configuration, the hinge is
% allowed to move resulting in a folding motion. 

% Initial x-values
x1 = [0; 0];
x2 = [1; 0];
x3 = [2; 0];

% x rigidity constraint matrix 
U = [-eye(2), zeros(2), eye(2)
    eye(2), zeros(2), zeros(2)
    zeros(2), eye(2), zeros(2)];

% constraining the origin and the length of the outer edges of the bar
A = [-eye(3), zeros(3), eye(3)
    eye(3), zeros(3), zeros(3)];

% final x vector
x = [x1; x2; x3];

% Initial y-values
y1 = [0; 0; 0];
y2 = [0.75; 1; 0];
y3 = [1.5; 0; 0];

% final y vector
y = [y1; y2; y3];

% populate the vector numbering all of the panels
J = [1, 2];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2];
F2 = [2, 3];

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

tol = 0.00001;

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "1DBar");


