function x_matrixSum = x_calcMatrixSum(zM, Fj)
% 
% calculates the summation of matrices Aij that allows the y vector to be
% factored out
%
% inputs:
% Fj: a cell array of the set of all y within a given panel 
% zM: cell array of matrices
%
%
% outputs:
% matrixSum: the sum of all of the matrices corresponding to a given panel
% 
% 
% testing to see if this changes anything
%
%

m = length(zM);
sum = zeros(3, 3*m);

for i = 1:length(Fj) 
    q = Fj{i};
    sum = sum + M{q};
end

x_matrixSum = x_sum; 

end