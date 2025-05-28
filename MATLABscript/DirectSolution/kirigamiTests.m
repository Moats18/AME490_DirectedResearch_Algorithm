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
% Updated Date: 03/13/25

% Initial x-values: planar kirigami, see Fig. S1 from "Programming bistability in
% geometrically perturbed mechanical metamaterials" Y. Peng et. al
%{
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
%}

% sheared planar kirigami, see Fig. 4 and Eq. (6) from above paper 
% gamma = 0.2; % Fig. 4 shows stability all the way up to gamma = 0.95
%l1D = l1R + gamma*l2R;
%l2D = l2R + gamma*l1R;
%{
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
%}

% Reference and deformed configuration lattice vector parameters
test = "shear";

s_r = 1; % reference side length
l1R = s_r*[1; 0];
l2R = s_r*[0; 1];

lambda_1 = 1.2; % axial deformation
lambda_2 = 0.8;
l1D_axial = lambda_1*l1R;
l2D_axial = lambda_2*l2R;

gamma = 0.90; % shear deformation; try 0.2, 0.35, 0.50, 0.65, 0.80, 0.95
l1D_shear = l1R + gamma*l2R;
l2D_shear = gamma*l1R + l2R;

if test == "reference"
    s = 0.5 * s_r;
elseif test == "axial"
    l1R = l1D_axial;
    l2R = l2D_axial;
    s_1 = lambda_1;
    s_2 = lambda_2;
elseif test == "shear"
    l1R = l1D_shear;
    l2R = l2D_shear;
    s = sqrt(s_r^2+gamma^2)/2;
end

% Initial x-values for planar kirigami initialized from lattice vectors & weighted averages 
x1 = [0; 0];
x5 = x1 + l1R;
x2 = 3/4 * x1 + 1/4 * x5;
x4 = 1/4 * x1 + 3/4 * x5;
x11 = x1 + l2R;
x12 = 1/2 * x11;
x7 = x5 + x11;
x6 = 1/2 * x5 + 1/2 * x7;
x10 = 3/4 * x11 + 1/4 * x7;
x8 = 1/4 * x11 + 3/4 * x7;
%x3 = 3/4 * 1/2 * (x6 + x12) + 1/4 * 1/2 * (x2 + x4);
%x9 = 3/4 * 1/2 * (x6 + x12) + 1/4 * 1/2 * (x8 + x10);
x3 = 3/4 * 1/2 * (x2 + x4) + 1/4 * 1/2 * (x10 + x8);
x9 = 3/4 * 1/2 * (x10 + x8) + 1/4 * 1/2 * (x2 + x4);

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

% Initial y-values for rotating squares with xi degree of separation
%xi = asin(x9(2) - x3(2));
xi = asin(0.9);

if test == "axial"
    y1 = [0; 0; 0];
    y2 = [s_1*cos(xi); s_2*-1*sin(xi); 0];
    y3 = [s_1sin(xi)+cos(xi); s_2*(cos(xi)-sin(xi)); 0];
    y4 = [s_1*2*sin(xi)+cos(xi); s_2*-1*sin(xi); 0];
    y5 = [s_1*2*sin(xi)+2*cos(xi); 0; 0];
    y6 = [s_1*sin(xi)+2*cos(xi); s_2*cos(xi); 0];
    y7 = [s_1*2*sin(xi)+2*cos(xi); s_2*2*cos(xi); 0];
    y8 = [s_1*2*sin(xi)+cos(xi); s_2*(2*cos(xi)+sin(xi)); 0];
    y9 = [s_1*sin(xi)+cos(xi); s_2*(cos(xi)+sin(xi)); 0];
    y10 = [s_1*cos(xi); s_2*(2*cos(xi)+sin(xi)); 0];
    y11 = [s_1*0; s_2*2*cos(xi); 0];
    y12 = [s_1*sin(xi); s_2*cos(xi); 0];
else
    y1 = [0; 0; 0];
    y2 = s*[cos(xi); -sin(xi); 0];
    y3 = s*[sin(xi)+cos(xi); cos(xi)-sin(xi); 0];
    y4 = s*[2*sin(xi)+cos(xi); -sin(xi); 0];
    y5 = s*[2*sin(xi)+2*cos(xi); 0; 0];
    y6 = s*[sin(xi)+2*cos(xi); cos(xi); 0];
    y7 = s*[2*sin(xi)+2*cos(xi); 2*cos(xi); 0];
    y8 = s*[2*sin(xi)+cos(xi); 2*cos(xi)+sin(xi); 0];
    y9 = s*[sin(xi)+cos(xi); cos(xi)+sin(xi); 0];
    y10 = s*[cos(xi); 2*cos(xi)+sin(xi); 0];
    y11 = s*[0; 2*cos(xi); 0];
    y12 = s*[sin(xi); cos(xi); 0];
end

y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12];

% x rigidity constraint matrix (5x12) of (2x2) = (10x24)

U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2)];
%{
U = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
%}
% y rigidity constraint matrix (5x12) of (3x3) = (15x36)

A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3)];
%{
A = [-eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];
%}
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
tol = 10^(-4);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");

