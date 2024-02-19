function pj = rotationMinimization(x, y, Fj, Tj)
%ROTATIONMINIMIZATION Based on fixed arrays of X and Y 
%
% inputs:
% takes the set of all x's within a panel (called Tj)
% all x coordinates 
% takes the set of all y's within a panel (called Fj)
% and all y coordinates 
%
% outputs:
% the quaternion pj that minimizes the energy with respect to a rotation 
%

size = length(Fj);
Bj = zeros(4, 4);

[cj, ~] = centerOfPanel(Fj, y);
[~, rij] = centerOfPanel(Tj, x);

for i = 1:size
V = rij(i) + cj - y{i};
T = cj -rij(i) - y{i};
Bij{i} = [0 V(1) V(2) V(3); T(1) 0 -V(3) V(2); T(2) V(3) 0 -V(1); T(3) -V(2) V(1) 0];
end 

for i = 1:size
Bj = Bj + Bij{i};
end 

[Vec ,~] = eig(Bj'*Bj);
pj = Vec(:, 1);

end

