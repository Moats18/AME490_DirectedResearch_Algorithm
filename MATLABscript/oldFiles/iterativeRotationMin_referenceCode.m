function Ropt = iterativeRotationMin(x, y, Fj, Tj, Ti, J, tol, R)
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
% Ti: a cell array containing the set of all panels associated with index i
% J: the set of all panels
% tol: the tolerance that the convergence has to be less than
% R: a cell array of the initial rotation matrix for each panel 
%
% outputs:
% Ropt: a cell array of the rotation matrices that minimize the elastic energy 

%initial values
E{1} = 0;
yNew = zeros(length(y), 1);
rij = zeros(3*length(Tj), 1, length(J));

% setting the values of the initial cj and rij vectors
for j = 1:length(J)
% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj(:, :, j), y);

% pos vectors with respect to the center of the panel
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

% calculating the elastic energy based on the initial cj, y, and x values
for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    k = Fj(:, i, j);
    E{1} = E{1} + norm(y(3*k-2:3*k, 1) - cj{j} - R{j}*rij(3*i-2:3*i, 1, j))^2;   
    end
end 

% setting the new values of the y position vector
for i = 1:(length(y)/3) 
    sum = zeros(3, 1);
    count = 1;
     j = Ti{i};
        for k = 1:length(Ti{i})
        l = j(k);
        [~, index] = find(Tj(:, :, l) == i);
        sum = sum + cj{l} + R{l}*rij(3*index-2:3*index, 1, l);% need to figure out R{j} vector
        count = count+1;
        end
    yNew(3*i-2:3*i, 1) = sum/length(Ti{i});
end


% updating the value of the center of the panel
for j = 1:length(J)
[cjNew{j}, ~] = centerOfPanel(Fj(:, :, j), yNew);
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
    minEig = min(D(D~=0));
    [~, col] = find(D == minEig);
    if length(col) > 1
        eigenVal = col(2);
    else
        eigenVal = col;
    end
    pj{j} = V(:, eigenVal);
    Rnew{j} = quat2rotm(pj{j}');
end

n = 2;
E{n} = 0; 

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    k = Fj(:, i, j);
    E{n} = E{n} + norm(yNew(3*k-2:3*k, 1)- cjNew{j} - Rnew{j}*rij(3*i-2:3*i, 1, j))^2;

    end
end 

R = Rnew; 
y = yNew;
cj = cjNew;

% need to update the new chnages
while (abs(E{n} - E{n-1}) > tol)

% setting the new values of the y position vector
for i = 1:(length(y)/3) 
    sum = 0;
    count = 1;
     j = Ti{i};
        for k = 1:length(Ti{i})
        l = j(k);
        [~, index] = find(Tj(:, :, l) == i);
        sum = sum + cj{l} + R{l}*rij(3*index-2:3*index, 1, l);% need to figure out R{j} vector
        count = count+1;
        end
    yNew(3*i-2:3*i, 1) = sum/length(Ti{i});
end

% updating the value of the center of the panel
for j = 1:length(J)
[cjNew{j}, ~] = centerOfPanel(Fj(:, :, j), yNew);
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
    minEig = min(D(D~=0));
    [~, col] = find(D == minEig);
    if length(col) > 1
        eigenVal = col(2);
    else
        eigenVal = col;
    end
    pj{j} = V(:, eigenVal);
    Rnew{j} = quat2rotm(pj{j}');
end

n = n+1;
E{n} = 0; 

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    k = Fj(:, i, j);
    E{n} = E{n} + norm(yNew(3*k-2:3*k, 1)- cjNew{j}- Rnew{j}*rij(3*i-2:3*i, 1, j))^2;
    end
end 

% Updating the values to be used at the beginning of the next loop
R = Rnew; 
y = yNew;
cj = cjNew;
end

Ropt = R;

end 