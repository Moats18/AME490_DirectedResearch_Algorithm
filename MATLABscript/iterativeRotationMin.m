function Ropt = iterativeRotationMin(x, y, Fj, Tj, Ti, J, tol, R)
% 
% uses quaternions to conver the problem into an eigenvalue probelm, 
% the iterative rotation minimization produces a rotation matrix that
% yields a convergent, minimized energy based on a configuration of
% constant X and Y values
%
% inputs:
% Tj: a cell array of the set of all x's within each panel 
%
% ***** update x structure***
% x: x coordinate 2-D array (3*n by 1 where n is the number of indices)
%
% Fj: a cell array of the set of all y's within each panel 

% ***** update y structure***
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
%
% Ti: the set of all panels associated with index i
% J: the set of all panels
% tol: the tolerance that the convergence has to be less than
% R: the rotation matrix for each panel 
%
% outputs:
% Ropt: rotation matrix that minimizes the elastic energy 

%initial values
Ri = eye(3);
E{1} = 0;

% setting the values of the initial cj and rij vectors
for j = 1:length(J)
% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj{j}, y);

% pos vectors with respect to the center of the panel
[~, rij{:, j}] = centerOfPanel(Tj{j}, x);
end 

% calculating the elastic energy based on the initial cj, y, and x values
for j = 1:length(J)
    for i = 1:length(Fj{j})
    E{1} = E{1} + norm(y{i} - cj{j} - Ri*rij{i, j})^2;
    end
end 

% setting the new values of the y position vector
for i = 1:length(y)   
    for k = 1:length(Ti)
        j = Ti{k};
        sum = cj{j} + R{j}*rij{i, j};% need to figure out R{j} vector
    end 
    yNew{i} = sum/length(Ti);
end

% updating the value of the center of the panel
[cjNew{j}, ~] = centerOfPanel(Fj{j}, yNew);

% creation of the rotation matrix Bij
for j = 1:length(J)
    for i = 1:length(y)
        % ***** update y structure***
       V = rij(i) + cj - y{i}; T = cj -rij(i) - y{i};
       Bij{i, j} = [0 V(1) V(2) V(3); T(1) 0 -V(3) V(2); T(2) V(3) 0 -V(1); T(3) -V(2) V(1) 0];
    end
end

for j = 1:length(J)
   size = length(Fj{j});
    for i = 1:size
    Bj{j} = Bj{j} + Bij{i, j};
    end 
    [Vec ,~] = eig(Bj{j}'*Bj{j});
    pj{j} = Vec(:, 1);
    Rnew{j} = quat2rotm(pj{j});
end

n = 2;
E{n} = 0; 

for j = 1:length(J)
    for i = 1:length(Fj{j})
    E{n} = E{n} + norm(yNew{i} - cjNew{j} - Rnew{j}*rij{i, j})^2;
    end
end 


while E{n} - E{n-1}>tol

% setting the new values of the y position vector
for i = 1:length(y)   
    for k = 1:length(Ti)
        j = Ti{k};
        sum = cj{j} + R{j}*rij{i, j};% need to figure out R{j} vector
    end 
    yNew{i} = sum/length(Ti);
end

% updating the value of the center of the panel
[cjNew{j}, ~] = centerOfPanel(Fj{j}, yNew);

% creation of the rotation matrix Bij
for j = 1:length(J)
    for i = 1:length(y)
       V = rij(i) + cj - y{i}; T = cj -rij(i) - y{i};
       Bij{i, j} = [0 V(1) V(2) V(3); T(1) 0 -V(3) V(2); T(2) V(3) 0 -V(1); T(3) -V(2) V(1) 0];
    end
end

for j = 1:length(J)
   size = length(Fj{j});
    for i = 1:size
    Bj{j} = Bj{j} + Bij{i, j};
    end 
    [Vec ,~] = eig(Bj{j}'*Bj{j});
    pj{j} = Vec(:, 1);
    Rnew{j} = quat2rotm(pj{j});
end

n = n+1;
E{n} = 0; 

for j = 1:length(J)
    for i = 1:length(Fj{j})
    E{n} = E{n} + norm(yNew{i} - cjNew{j} - Rnew{j}*rij{i, j})^2;
    end
end 

% Updating the values to be used at the beginning of the next loop
R = Rnew; 
y = yNew;
cj = cjNew;
end

Ropt = R;

end 
