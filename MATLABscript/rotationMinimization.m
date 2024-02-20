function pj = rotationMinimization(x, y, Fj, Tj)
%ROTATIONMINIMIZATION Based on fixed arrays of X and Y 
%
% inputs:
% Tj: takes the set of all x's within a panel
% x: 2-D array of all of the x coordinates (3*n by 1 where n is the number of indices)
% Fj: takes the set of all y's within a panel 
% y: 2-D array of all of the y coordinates (3*n by 1 where n is the number of indices)
%
% outputs:
% pj: a column vector that includes the quaternion pj that minimizes the energy with respect to a rotation 
%

size = length(Fj);
Bj = zeros(4, 4);

[cj, ~] = centerOfPanel(Fj, y);
[~, rij] = centerOfPanel(Tj, x);

for i = 1:size
V = rij(3*i-2:3*i, 1) + cj - y(3*i-2:3*i, 1);
T = cj - rij(3*i-2:3*i, 1) - y(3*i-2:3*i, 1);
Bij{i} = [0 V(1) V(2) V(3); T(1) 0 -V(3) V(2); T(2) V(3) 0 -V(1); T(3) -V(2) V(1) 0];
end 

for i = 1:size
Bj = Bj + Bij{i};
end 

[Vec ,~] = eig(Bj'*Bj);
pj = Vec(:, 1);

end

