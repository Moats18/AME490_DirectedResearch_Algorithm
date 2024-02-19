posVec = [0,0,0, 1,0,0, 0,1,0, 1,1,0, 2,1,0, 2,2,0]';
numPanels = 2;
dim = 3;
len = length(posVec);

panel{1} = [1, 2, 3, 4];
panel{2} = [3, 4, 5, 6];

% create the x matrices
for k = 1:6
    x{k} =  [zeros(dim, 3*(k-1)), eye(dim), zeros(dim, 3*(6-k))];
   
end

lenX = length(x);

% creating the center of panel matrix 
for j = 1:numPanels
    sum = zeros(3, 18);
    for k = panel{j}
    sum = sum + x{k}; % summing over all the matrices given a panel 
    end
    c{j} = 1/4*sum; 
end

k = zeros(18,18);

%k matrix
for j = 1:numPanels
    sum = zeros(3, 18);
    for i = panel{j}
        sum = sum + x{i} - c{j}; % equal to Aij
        k = k + sum'*sum; % equal to Aij transpose * Aij
    end
end


% b matrix 

