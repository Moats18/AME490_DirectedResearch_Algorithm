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
% Updated Date: 05/23/24
%
% The initial configuration consists of four panels that share common
% edges. All of the nine vertices of the panels can be determined based on two user 
% specified points. Using the symmetry constraints imposed on these two points 
% results in nine fully defined points

% Initial x-values
x1 = [0; 0; 0];
x2 = [1; 1; 0];

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
    
% setting the x and y constraint matrix to be the same i.e. stating that
% they have the same overall shape
U = A;

% x lattice vectrs
e1 = [0; 2;0];
e2 = [1;0;0];
e3 = [0;0;0]; % new addition
h = [e2; e1; e2; e2; e3; e2; e2; e2]; 

% populate the x vector
xU = A(1:24, 7:27)\(h-(A(1:24, 1:6)*[x1; x2])); % 

% final x vector
x = zeros(27, 1); % chnaged size from 18 to 27
x(1:6, 1) = [x1; x2];
x(7:27, 1) = xU; 

% Initial y-values
y1 = [0;0;0];
y2 = [0; 1;0];

% y lattice vectors
l1 = [0; 2; 0];
l2 = [1; 0; 0];
l3 = [0;0;0];% new addition
e = [l2; l1; l2; l2; l3; l2; l2; l2]; 

% populate the y vector
yU = A(1:24, 7:27)\(e-(A(1:24, 1:6)*[y1; y2]));% chnaged size from 1:15 to 1:18

% final y vector
y = zeros(27, 1);% chnaged size from 15 to 18
y(1:6, 1) = [y1; y2];
y(7:27, 1) = yU;% chnaged size from 1:15 to 1:18

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

% determine the index set for each index
t1 = 1;
t2 = [1, 2];
t3 = 2;
t4 = [2, 3];
t5 = [1, 2, 3, 4];
t6 = [1, 4];
t7 = 4;
t8 = [3, 4];
t9 = 3;
Ti = {t1, t2, t3, t4, t5, t6, t7, t8, t9}; %the set of all panels associated with index i

% Initial R 
for j = 1:length(J)
R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for the entire algorithm minimization
tol = 0.001;

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = [x, y, xOpt, yOpt];
visualizeLatticeVec = true;

figure
plot4vectors(vectors, titles, visualizeLatticeVec);

