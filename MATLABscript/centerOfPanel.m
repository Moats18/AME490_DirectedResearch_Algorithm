function [cj, rij] = centerOfPanel(Fc, coord)
% Inputs:
% Fc: The set of all coordinates (x or y) within a panel 
%
% ***** update coord structure***
% coord: 2-D array of all of the x/y coordinates (3*n by 1 where n is the number of indices)
%
% Outputs:
% cj: the center of the panel 
% rij: the coordinate vectors with respect to the center of the panel

len = length(Fc);
rij = zeros(3, len);
sum = zeros(1, 3);
num = 0;

for i = 1:length(Fj)
sum = sum + coord{Fc(i)};
num = num + 1;
end

cj = sum/num;

for i = 1:len
rij(:, i) = y{i} - cj;
end

end

