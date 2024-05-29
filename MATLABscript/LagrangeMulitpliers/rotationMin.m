function Ropt = rotationMin(x, y, Fj, Tj, J)
% 
% uses quaternions to conver the minimization into an eigenvalue probelm, 
% the iterative rotation minimization determines the rotation matrix for each panel
% that yields a convergent, minimized energy based on a configuration of
% given X and Y values
%
% inputs:
% Tj: a cell array of the set of all x's within each panel 
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
% Fj: a cell array of the set of all y's within each panel 
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
%
% J: the set of all panels
% tol: the tolerance that the convergence has to be less than
% R: a cell array of the initial rotation matrix for each panel 
%
% outputs:
% Ropt: a cell array of the rotation matrices that minimize the elastic energy 

% setting the values of the initial cj and rij vectors
for j = 1:length(J)

% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj(:, :, j), y);

% pos vectors with respect to the center of the panel
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

% creation of the rotation matrix Bij
for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
       k = Fj(1, i, j); 
       V = rij(3*i-2:3*i, 1, j) + cj{j} - y(3*k-2:3*k, 1);
       T = cj{j} - rij(3*i-2:3*i, 1, j) - y(3*k-2:3*k, 1);
       Bij{i, j} = [0 V(1) V(2) V(3); T(1) 0 -V(3) V(2); T(2) V(3) 0 -V(1); T(3) -V(2) V(1) 0];
    end
end

for j = 1:length(J)
Bj{j} = zeros(4);
end

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    Bj{j} = Bj{j} + Bij{i, j}'*Bij{i, j};
    end

    % rounding the small numeric values to zero
    roundedBj = Bj{j};
    roundedBj(abs(roundedBj)<1e-3)=0;

    [V ,D] = eig(roundedBj);
    roundedD = D;
    roundedD(abs(roundedD)<1e-3)=0;
    maxEig = max(roundedD, [], "all");
    [~, col] = find(roundedD == maxEig);

    if length(col) > 1
        eigenVal = col(1);
    else
        eigenVal = col;
    end

    pj{j} = V(:, eigenVal);
    Rnew{j} = quat2rotm(pj{j}');
end

Ropt = Rnew;

end 