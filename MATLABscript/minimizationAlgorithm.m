function [yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, J, R, A, e. tol)
% 
% Based on the paper: "Elastic Energy Approximation and Minimization
% Algorithm for Foldable Meshes" 
%
% By: Antoine Moats, Niharika Sashidhar
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Date: 02/24/24
%
% Rigidity Constraints:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
% e: e is a vector of the rigid lengths between different y coordinates
%
% Indexing Inputs:
% Tj: a cell array of the set of all x's within each panel(the jth panel
% corresponds to the jth row)
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a 2D array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% J: the set of all panels
% R: a cell array of all of the initial rotation matrix for each panel
%
% Outputs:
% yMin: y coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% xMin: x coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% Ropt: array of rotation matrices for minimizes the elastic energy  

% initial values
err = inf;
E{1} = 0;

% setting the values of the initial cj and rij vectors
for j = 1:length(J)
% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj{j}, y);

% pos vectors with respect to the center of the panel
[~, rij{:, :, j}] = centerOfPanel(Tj{j}, x);
end 

for j = 1:length(J)
    for i = 1:length(Fj{j})
    E{1} = E{1} + norm(y{i} - cj{j} - Ri*rij{:, i, j})^2;
    end
end 

while err > tol 



err = abs(E{n}-E{n-1});

end




end

