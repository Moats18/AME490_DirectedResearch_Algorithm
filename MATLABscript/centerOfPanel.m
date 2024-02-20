function [cj, rij] = centerOfPanel(Fc, coord)
% Inputs:
% Fc: The set of all coordinates (x or y) within a panel 
% coord: 2-D array of all of the x/y coordinates (3*n by 1 where n is the number of indices)
%
% Outputs:
% cj: the center of the panel 
% rij: the coordinate vectors with respect to the center of the panel

len = length(Fc);
rij = zeros(3*len, 1);
sum = zeros(3, 1);
num = 0;

for i = 1:length(Fc)
sum = sum + coord(3*Fc(i)-2:3*Fc(i), 1);
num = num + 1;
end

cj = sum/num;

for i = 1:len
rij(3*i-2:3*i, 1) = coord(3*Fc(i)-2:3*Fc(i), 1) - cj;
end

end

