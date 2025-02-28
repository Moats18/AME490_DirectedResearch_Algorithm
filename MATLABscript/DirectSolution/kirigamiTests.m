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
% Updated Date: 02/28/25

% Initial x-values: planar kirigami, see Fig. S1 from "Programming bistability in
% geometrically perturbed mechanical metamaterials" Y. Peng et. al
x1 = [0; 0];
x2 = [0.5; -1];
x3 = [0.4; 3];
x4 = [1.3; -1.1];
x5 = [2.6; 1.5];
x6 = [2; 4];
x7 = [2.2; 8];
x8 = [1.6; 7.5];
x9 = [0.7; 5];
x10 = x8+(x2-x4); % constrained for tessellation
x11 = x1+(x7-x5); % constrained for tessellation
x12 = [-0.3; 6];

l1R = x5-x1; % lattice vectors
l2R = x10-x2;

% sheared planar kirigami, see Fig. 4 and Eq. (6) from above paper
gamma = 0.2; % Fig. 4 shows stability all the way up to gamma = 0.95
l1D = l1R + gamma*l2R;
l2D = l2R + gamma*l1R;

x1 = [0; 0];
x2 = [1; 0];
x3 = [1; 1];
x4 = [1.5; 0.5];
x5 = [2; 1];
x6 = [1.8; 1.5];
x7 = [2.5; 2.5];
x8 = [2.2; 3];
x9 = [1.25; 1.25];
x10 = x8+(x2-x4); % constrained for tessellation
x11 = x1+(x7-x5); % constrained for tessellation
x12 = [0.3; 0.5];

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

% Initial y-values
s = 1;
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
y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12];

% x rigidity constraint matrix (4x12) of (2x2) = (8x24)
U = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2)];

% y rigidity constraint matrix (4x12) of (3x3) = (12x36)
A = [-eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3);
     zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];

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

% 3D array that contains the index set of x coordinates within a given
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
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");

