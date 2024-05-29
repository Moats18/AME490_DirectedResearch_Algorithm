function E = calcEnergy(vec, Tj, Fj, J)

x = vec(1:27);
y = vec(28:54);

E = 0;
rij = zeros(3*length(Tj), 1, length(J));

% setting the values of the initial cj and rij vectors

for j = 1:length(J)
% center of the panel calculation based on initial y vector
[cj{j}, ~] = centerOfPanel(Fj(:, :, j), y);

% pos vectors with respect to the center of the panel
[~, rij(:, :, j)] = centerOfPanel(Tj(:, :, j), x);
end 

for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
    k = Fj(:, i, j);
    E = E + norm(y(1, 3*k-2:3*k) - cj{j} - rij(3*i-2:3*i, 1,j)')^2; 
    end
end 


end

