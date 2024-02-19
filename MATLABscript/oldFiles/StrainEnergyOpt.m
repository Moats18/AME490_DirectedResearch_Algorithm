
%initializing the position vectors 

% for a mesh of maximum y value and maximum x value 
maxY = 2;
maxX = 1;
index = 1;
numPanels = maxY*maxX;
dim = 3;

for i = 0:maxY 
    for j = 0:maxX
    posVec{index} = [j, i, 0];
    index = index + 1;
    end
end

len = length(posVec);

% creating the panel indexes
for p = 1:numPanels
        for l 
        panelInd{p, l} = ;
end

% 


% create the x matrices
for k = 1:len
    x{k} =  [zeros(dim, k-1), eye(dim), zeros(dim, len-k)];
end

%create the A matrix






