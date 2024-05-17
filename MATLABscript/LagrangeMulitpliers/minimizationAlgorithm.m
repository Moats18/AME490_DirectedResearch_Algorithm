function [yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, J, R, A, U, h, e, tol)
% 
% Based on the paper: "Elastic Energy Approximation and Minimization
% Algorithm for Foldable Meshes" 
%
% By: Antoine Moats, Niharika Sashidhar
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Date: 03/27/24
%
% Rigidity Constraints:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
% e: e is a vector of the rigid lengths between different y coordinates
% U: U is a matrix that contains the rigidity constraints such that Ux = h
% h: h is a vector of the rigid lenghts between different x coordinates\
%
% Indexing Inputs:
% Tj: a 3D array of the set of all x's within each panel 
% Tj size is 1 by verticesPerPanel  by nuPanels;
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a 2D array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% J: the set of all panels
% R: a cell array of all of the initial rotation matrix for each panel
% tol: tolerance for the algorithm minimization
%
% Outputs:
% yMin: y coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% xMin: x coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% Ropt: array of rotation matrices for minimizes the elastic energy  

% initial values
loop = 1; % loop number
E{1} = 0;
rij = zeros(3*length(Tj), 1, length(J));

% setting the values of the initial cj and rij vectors
for j = 1:length(J)
% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj(:, :, j), y);

% pos vectors with respect to the center of the panel
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    k = Fj(:, i, j);
    E{1} = E{1} + norm(y(3*k-2:3*k, 1) - cj{j} - R{j}*rij(3*i-2:3*i, 1, j))^2;
    end
end 

% defining the first index that is to be compared in the while loop
% conditional statement
num = 2;
E{num} = 0; 

% Calculating the first energy value
n = length(y);
lenJ = length(J);

for i = 1:n/3
    xX = zeros(3, n);
    xX(1:3, 3*i-2:3*i) = eye(3);
    xM{i} = xX;
end 

% the first minimized rotation matrix
RiOpt = iterativeRotationMin(x, y, Fj, Tj, J);

% the first minimized y 
yNew = minY(x, y, Fj, Tj, J, RiOpt, A, e);

figure
subplot3dvec(y, 'y', yNew, 'yNew');

% the first minimized x 
xNew = minX(x, yNew, Fj, Tj, J, RiOpt, U, h);

figure
subplot3dvec(x, 'x', xNew, 'xNew')
hold off 

% Defining the matrix that allows the vector y to be factored out
% Aij{1, 1} is a 3 by 3*n matrix
for j = 1:lenJ
    l = length(Fj(:, :, j));
    sum = (1/l)*calcMatrixSum(xM, Fj(:, :, j));
    for i = 1:length(Fj(:, :, j))
        k = Fj(:, i, j);
        Aij{i, j} = xM{k} - sum;
    end
end

% Defining the matrix that allows the vector x to be factored out
% Qij{1, 1} is a 3 by 3*n matrix
for j = 1:lenJ
    l = length(Fj(:, :, j));
    sum = (1/l)*calcMatrixSum(xM, Fj(:, :, j));
    for i = 1:length(Fj(:, :, j))
        k = Fj(:, i, j);
        Qij{i, j} = xM{k} - sum;
    end
end

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    E{num} = E{num} + norm(Aij{i, j}*yNew - RiOpt{j}*(Qij{i, j}*xNew))^2;
    end
end 

R = RiOpt; 
y = yNew;
x = xNew;

err = abs(E{num}-E{num-1});

% while loop that converges on a minimized energy value
figure
xlabel("Iteration Number [#]")
ylabel("Energy [m^2]");
title("Minimization Algorithm");

hold on
for i=2:length(E)
scatter(i, E{i}, 'filled', 'MarkerFaceColor', [0.10, 0.60, 0.9]);
drawnow; % ensures that the updated point is plotted
pause(0.5); %pausing for 1/10 of a second
end

while err > tol 

Ropt = iterativeRotationMin(x, y, Fj, Tj, J);
yNew = minY(x, y, Fj, Tj, J, Ropt, A, e); 
xNew = minX(x, yNew, Fj, Tj, J, Ropt, U, h);

count = num +  loop;
E{count} = 0;
for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    E{count} = E{count} + norm(Aij{i, j}*yNew - Ropt{j}*(Qij{i, j}*xNew))^2;
    end
end 

err = abs(E{count}-E{count-1});

% Updating the values to be used at the beginning of the next loop
R = Ropt; 
y = yNew;
x = xNew;

% showing the results of the algorithm in real time
disp("-------------------------------------")
disp("Iteration number: " + num2str(loop));
disp("Total Energy Calculations: " + num2str(count));
disp("-------------------------------------")
disp("Previous Energy Value: " + num2str(E{count-1}));
disp("Current Energy Value: " + num2str(E{count}));
disp("Current Error Value: " + num2str(err));

% plotting the results of the algorithm in real time
hold on
scatter(count, E{count}, 'filled', 'MarkerFaceColor', [0.10, 0.60, 0.9]);
plot(count, E{count});
drawnow; % ensures that the updated point is plotted
pause(0.5); %pausing for 1/2 of a second
loop = loop + 1;
hold off
end

Ropt = R;
yOpt = y;
xOpt = x;

end

