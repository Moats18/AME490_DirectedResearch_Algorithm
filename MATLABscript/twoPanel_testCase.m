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
% Date: 03/27/24
%
% the initial configuration consists of two panels that share a common edge
% all of the six vertices of the panel can be determined based on two randomly 
% generated points. Using the symmetry constraints imposed on these two points 
% results in six fully defined points

% Initial x-values
x1 = randi(10, 3, 1); %generates an array of 3 random numbers from 1-10 
x2 = randi(10, 3, 1); %generates an array of 3 random numbers from 1-10 

% x rigidity constraints- define as 
A = [-eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),eye(3,3);
    -eye(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3);
    zeros(3,3),-eye(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3);
    zeros(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3), -eye(3, 3);
    zeros(3,3),zeros(3,3),-eye(3,3),eye(3,3),zeros(3,3),zeros(3,3)];

%
U = A;

% x - rigidity constraints
e1 = randi(5, 3, 1);
e2 = randi(5, 3, 1);
h = [e1; e2; e1; e2; e1];

% populate the x vector
xU = A(1:15, 7:18)\(h-(A(1:15, 1:6)*[x1; x2]));

% final x vector
x = zeros(15, 1);
x(1:6, 1) = [x1; x2];
x(7:18, 1) = xU;

% y - rigidity constraints
l1 = randi(5, 3, 1);
l2 = randi(5, 3, 1);
e = [l1; l2; l1; l2; l1];

% Initial y-valudes
y1 = randi(10, 3, 1); %generates an array of 3 random numbers from 1-10 
y2 = randi(10, 3, 1); %generates an array of 3 random numbers from 1-10 

% populate the y vector
yU = A(1:15, 7:18)\(e-(A(1:15, 1:6)*[y1; y2]));

% final y vector
y = zeros(15, 1);
y(1:6, 1) = [y1; y2];
y(7:18, 1) = yU;

% populate the vector numbering all of the panels
J = [1, 2];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
% a cell array of the set of all x's within each panel 
T1 = [1, 2, 5, 6];
T2 = [2, 3, 4, 5];

% a 3-d array that contains the index set of x coordinates within a given
% panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;

Fj = Tj;

% determine the index set for each index
t1 = 1;
t2 = [1, 2];
t3 = 2;
t4 = 2;
t5 = [1, 2];
t6 = 1;
Ti = {t1, t2, t3, t4, t5, t6}; %the set of all panels associated with index i

% Initial R 
for j = 1:length(J)
R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for the rotation minimization
tolR = 0.1;

% determine the initial tolerance for the entire algorithm minimization
tol = 0.1;

[yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, Ti, J, R, A, U, h, e, tol, tolR);

