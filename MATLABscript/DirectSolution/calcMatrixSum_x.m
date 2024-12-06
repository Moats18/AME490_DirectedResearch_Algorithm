function matrixSum = calcMatrixSum_x(xM, Tj)
% 
% calculates the summation of matrices Aij that allows the x vector to be
% factored out
%
% inputs:
% Tj: a cell array of the set of all x within a given panel 
% xM: cell array of matrices
%
%
% outputs:
% matrixSum: the sum of all of the matrices corresponding to a given panel
%

n = length(xM);
sum = zeros(3, 2*n);

for k = 1:length(Tj)
   p = Tj(k);
   sum = sum + xM{p};
end

matrixSum = sum; 

end

