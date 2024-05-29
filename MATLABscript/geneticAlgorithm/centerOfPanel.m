function [cj, rij] = centerOfPanel(Fc, coord)
% Inputs:
% Fc: The set of all of the index (as opposed to position) of coordinates (x or y) within a panel
% 1 by 4 2D array
% coord: 2D array of all of the x/y coordinates (3*n by 1 where n is the number of indices)
%
% Outputs:
% cj: the center of the panel 
% rij: 3 dimensional array containing the coordinate vectors with respect to the center of the panel
% ordered sequentially so iterating through i works 

len = length(Fc);
rij = zeros(1, 3*len);
sum = zeros(1, 3);
num = 0;

for i = 1:length(Fc)

sum = sum + coord(1, 3*Fc(i)-2:3*Fc(i));
num = num + 1;
end

cj = sum/num;

for i = 1:len
rij(1, 3*i-2:3*i) = coord(1, 3*Fc(i)-2:3*Fc(i)) - cj;
end

end

