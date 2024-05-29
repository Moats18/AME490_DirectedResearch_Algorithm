%
% This script is designed as a comparison of the built-in MATLAB genetic
% algorithm and the algorithm based on the paper:
% 
% % "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Antoine Moats, Niharika Sashidhar
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 05/23/24
%
% the initial configuration consists of two panels that share a common edge
% all of the six vertices of the panel can be determined based on two randomly 
% generated points. Using the symmetry constraints imposed on these two points 
% results in six fully defined points
%

% x rigidity constraint matrix
A = [-eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
    -eye(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
    zeros(3,3),-eye(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
    %zeros(3,3),zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3), -eye(3, 3); %c
    zeros(3,3),zeros(3,3),-eye(3,3),eye(3,3),zeros(3,3),zeros(3,3) ,zeros(3,3),zeros(3,3),zeros(3,3)
    eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
    zeros(3,3),zeros(3,3),zeros(3,3), -eye(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3), eye(3,3) %
    zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),-eye(3,3),zeros(3,3),zeros(3,3), eye(3,3), zeros(3,3) % 
    zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),-eye(3,3),eye(3,3),zeros(3,3), zeros(3,3)]; 

U = A;

% x lattice vectors
e1 = [0; 1.65; 0];
e2 = [0.75;0;0];
e3 = [0;0;0]; % new addition
h = [e2; e1; e2; e2; e3; e2; e2; e2]; 

% y lattice vectors
l1 = [0; 2; 0];
l2 = [1; 0; 0];
l3 = [0;0;0];% new addition
e = [l2; l1; l2; l2; l3; l2; l2; l2]; 

bigMatrix = zeros(48, 54);
bigMatrix(1:24, 1:27) = A;
bigMatrix(25:48, 28:54) = U;

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
F3 = [5, 4, 9, 8];
F4 = [6, 5, 8, 7];
% a cell array of the set of all x's within each panel 
T1 = [1, 2, 5, 6];
T2 = [2, 3, 4, 5];
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

% Initial x-values
x1 = [0; 0; 0];
x2 = [1; 0.75; 0];

% Initial y-values
y1 = [0;0;0];
y2 = [0; 1;0];

J = [1, 2, 3, 4];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
F3 = [5, 4, 9, 8];
F4 = [6, 5, 8, 7];
% a cell array of the set of all x's within each panel 
T1 = [1, 2, 5, 6];
T2 = [2, 3, 4, 5];
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

%genetic algorithm conditions

init = initialCond(x1, x2, y1, y2, A, U, e, h);
options = optimoptions('ga', 'InitialPopulationMatrix', init', 'PlotFcn', @gaplotbestf, 'MaxGenerations', 50);
fun = @(vec) calcEnergy(vec, Tj, Fj, J);
[optimizedVec, E] = ga(fun, 54, [], [], bigMatrix, [h, e], [], [], [], options);

%outputs
xopt = optimizedVec(1:27);
yopt = optimizedVec(28:54);
x = init(1:27)';
y = init(28:54)';
titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = [x, y, xopt, yopt];
visualizeLatticeVec = true;
plot4vectors(vectors, titles, visualizeLatticeVec);
disp("Final Energy: " + num2str(E));


