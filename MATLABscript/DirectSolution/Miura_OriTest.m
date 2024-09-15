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
% Updated Date: 08/12/24
%
% The initial configuration consists of four panels that share common
% edges. All of the nine vertices of the panels can be determined based on two user 
% specified points. Using the symmetry constraints imposed on these two points 
% results in nine fully defined points

% Initial x-values
s = 0.5;
gamma = pi/6;
x1 = [0; 0; 0];
x2 = [s*cos(gamma); s*sin(gamma); 0];
x3 = [2*s*cos(gamma); 0; 0];
x4 = [2*s*cos(gamma); s; 0];
x5 = [s*cos(gamma); s*sin(gamma)+s; 0];
x6 = [0; s; 0];
x7 = [0; 2*s; 0];
x8 = [s*cos(gamma); s*sin(gamma)+2*s; 0];
x9 = [2*s*cos(gamma); 2*s; 0];

% x rigidity constraint matrix (4x12) of (3x3) = (12x36)
U = [-eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     -eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3), -eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), zeros(3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3), eye(3,3)];
    
    
% setting the x and y constraint matrix to be the same i.e. stating that
% they have the same overall shape
A = [-eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
     -eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3);
     zeros(3,3), zeros(3), -eye(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), eye(3,3);
     zeros(3,3), zeros(3,3), zeros(3,3), zeros(3), zeros(3,3), zeros(3,3), -eye(3), zeros(3,3), eye(3,3)];

% final x vector
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

% Initial y-values
s = 0.5;
theta = pi/4;

% useful parameters of a Miura-Ori Origami
H = s*sin(theta)*sin(gamma);
S = s*cos(theta)*tan(gamma)/sqrt(1+cos(theta)^2*tan(gamma)^2);
L = s*sqrt(1-sin(theta)^2*sin(gamma)^2);
V = s*1/sqrt(1+cos(theta)^2*tan(gamma)^2);
eta = atan(cos(theta)*tan(gamma));
psi = asin(sin(theta)*sin(gamma));
phi = asin(sin(eta)/sin(gamma));
o = H/tan(theta);

y1 = lambda*[0; 0; 0];
y2 = lambda*[S; s*cos(eta); 0];
y3 = lambda*[2*S; 0; 0];
y4 = lambda*[2*S; L; H];
y5 = lambda*[S; 2*L; H];
y6 = lambda*[0; L; H];
y7 = lambda*[0; 2*L; 0];
y8 = lambda*[S; 2*L+V; 0];
y9 = lambda*[2*S; 2*L; 0];

% final y vector
y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];

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
tol = 10^(-2);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = [x, y, xOpt, yOpt];
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec);
