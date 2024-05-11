function yMin = minY_V2(x, y, Fj, Tj, J, R, A, e)
% 
% directly solves the linear equation to calculate the perturbation of the
% initial y 
%
% Inputs
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
% R: a cell array of all of the rotation matrices for each panel
%
% Outputs
% yMin: y coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   
% Define the matrices that allow the 3D coord vector yi to be converted to
% a list of coord vectors called y




end

