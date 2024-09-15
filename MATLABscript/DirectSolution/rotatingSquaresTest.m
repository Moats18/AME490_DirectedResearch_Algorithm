%
% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 08/12/24
%
% The initial configuration consists of four panels that share common
% edges. This constraint is relaxed following the rotating squares problem.
%
% Initial x-values
s = 0.5;
x1 = [0; 0; 0];
x2 = [s; 0; 0];
x3 = [s; s; 0];
x4 = [s; 0; 0];
x5 = [2*s; 0; 0];
x6 = [2*s; s; 0];
x7 = [2*s; 2*s; 0];
x8 = [s; 2*s; 0];
x9 = [s; s; 0];
x10 = [s; 2*s; 0];
x11 = [0; 2*s; 0];
x12 = [0; s; 0];

% x rigidity constraint matrix (4x12) of (3x3) = (12x36)
U = [-eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3);
     zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];
    
    
% setting the x and y constraint matrix to be the same i.e. stating that
% they have the same overall shape
A = [-eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3);
     zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];

% x lattice vectors
e1 = [0; 2*s; 0];
e2 = [2*s; 0; 0];
e3 = [0;0;0]; % new addition
h = [e2; e1; e2; e2; e3; e2; e2; e2]; 

% final x vector
x = zeros(3*12, 1);
x(1:3*12, 1) = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

% Initial y-values
lambda = sqrt(2);
y1 = lambda*[0; 0; 0];
y2 = lambda*[s; 0; 0];
y3 = lambda*[s; s; 0];
y4 = lambda*[s; 0; 0];
y5 = lambda*[2*s; 0; 0];
y6 = lambda*[2*s; s; 0];
y7 = lambda*[2*s; 2*s; 0];
y8 = lambda*[s; 2*s; 0];
y9 = lambda*[s; s; 0];
y10 = lambda*[s; 2*s; 0];
y11 = lambda*[0; 2*s; 0];
y12 = lambda*[0; s; 0];

% y lattice vectors
l1 = lambda*e1;
l2 = lambda*e2;
l3 = [0;0;0];% new addition
e = [l2; l1; l2; l2; l3; l2; l2; l2];

% final y vector
y = zeros(3*12, 1);
y(1:3*12, 1) = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12];

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 3, 12];
F2 = [3, 4, 5, 6];
F3 = [6, 7, 8, 9];
F4 = [9, 10, 11, 12];
% a cell array of the set of all x's within each panel 
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;

% a 3-d array that contains the index set of x coordinates within a given
% panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Tj(1, :, 3) = T3;
Tj(1, :, 4) = T4;

Fj = Tj;

% Initial R 
for j = 1:length(J)
R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = [x, y, xOpt, yOpt];
visualizeLatticeVec = true;

figure
plot4vectors(vectors, titles, visualizeLatticeVec);

