%
% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Antoine Moats, Niharika Sashidhar
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Date: 02/24/24
%

% Initial x-values

x{1} = [];
x{6} = [];
for i=2:5
    x{i} = randi(100,3,1);
end

% x rigidity constraints- define as 
A = {-eye(3,3),eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3);
    -eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3);
    zeros(3,3),-eye(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3);
    zeros(3,3),zeros(3,3),zeros(3,3),-eye(3,3),eye(3,3),zeros(3,3);
    zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3),-eye(3,3)};
    
% Initial y-valudes
y1 = [];
y6 = [];

% y rigidity constraints

% populate the x vector

% populate the y vector

% populate the vector numbering all of the panels
J = [];

% determine the index set for each panel 
Fj = {}; % a cell array of the set of all y's within each panel 
Tj = {}; % a cell array of the set of all x's within each panel 

% determine the index set for each index
Ti = []; %the set of all panels associated with index i

% Initial R 
R = eye(3); % identity matrix

% determine the initial tolerance for the rotation minimization
tolR = 0.1;

% determine the initial tolerance for the entire algorithm minimization
tol = 0.1;

[yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, J, R, A, U, h, e, tol, tolR);

