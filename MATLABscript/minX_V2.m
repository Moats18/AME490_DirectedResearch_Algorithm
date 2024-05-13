function Xmin = minX_V2(x, y, Fj, Sj, J, R, U,h)
%
% directly solves the linear equation to calculate the perturbation of the
% initial x
%
%Inputs:
%
% Rigidity Constraints:
% U: U is a matrix that contains the rigidity constraints such that Ux = h
% h: h is a vector of the rigid lenghts between different x coordinates\
%
% Indexing Inputs:
% x: x coordinate 2-D array (3*m by 1 where m is the number of indices)
% y: y coordinate 2-D array (3*m by 1 where m is the number of indices)
% Fj: a cell array of the set of all y's wihtin each panel
% Sj: a cell array of the set of all y's within each panel
% J: the set of all panels
% R: a cell array of all of the rotation matrices for each panel
% 
% Outputs:
%Xmin: x coordinate(2-D array) that minimizes the energy based on given
%constraints

% Define the matrices that allow the singular vector xi to be converted to
% a list of vectors called x


%creating the g matrix 


%


end

