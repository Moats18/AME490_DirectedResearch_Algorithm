function matrixSum = calcMatrixSum(xM, Fj)
% 
% calculates the summation of matrices Aij that allows the y vector to be
% factored out
%
% inputs:
% Fj: a cell array of the set of all y within a given panel 
% xM: cell array of matrices
%
%
% outputs:
% matrixSum: the sum of all of the matrices corresponding to a given panel
% 
% 
% testing to see if this changes anything
%
%

n = length(xM);
sum = zeros(3, 3*n);

for k = 1:length(Fj) 
    p = Fj{k};
    sum = sum + xM{p};
end

matrixSum = sum; 

end

