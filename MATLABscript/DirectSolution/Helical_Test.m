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
% Updated Date: 12/06/24.
%
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.

% Initial x-values
x1 = [0; 0];
x2 = [1; 0];
x3 = [2; 0];
x4 = [2; 1];
x5 = [1; 1];
x6 = [0; 1];
x7 = [0; 2];
x8 = [1; 2];
x9 = [2; 2];
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

p = 6;
q = 4;

tau1 = 0.5;
e1 = [0,0,1]';
z = [0,1,0]';
theta1 = pi/4;

[T1, T2, R1, R2] = calc_heli(q, p, theta1, tau1, z, e1);


% Initial y-values
y1 = [0; 0; 0];
y3 = R1*y1 + T1;
y7 = R2*y1 + T2;
y9 = R2*y3 + T2;

% average of corners
y2 = (y1+y3)/2;
y4 = (y3+y9)/2;
y6 = (y1+y7)/2;
y8 = (y7+y9)/2;

%average of middle points
y5 = (y6+y4)/2;

y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];

% symmetry constraints for reference
R1_x = eye(2);
R2_x = eye(2);
T1_x = zeros(2, 1);
T2_x = zeros(2, 1);

% determine constraint matrix U
U = convert_g2A2D(R1_x, R2_x, T1_x, T2_x);

% symmetry constraints for deformed
R1_y = eye(3);
R2_y = eye(3);
T1_y = zeros(3, 1);
T2_y = zeros(3, 1);

% determine constraint matrix U
A = convert_g2A3D(R1_y, R2_y, T1_y, T2_y);

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
F3 = [5, 6, 7, 8];
F4 = [4, 5, 8, 9];
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
tol = 10^(-3);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Helical");
