function init = initialCond(x1, x2, y1, y2, A, U, e, h)

% populate the y vector
yU = A(1:24, 7:27)\(e-(A(1:24, 1:6)*[y1; y2]));

% populate the x vector
xU = U(1:24, 7:27)\(h-(U(1:24, 1:6)*[x1; x2])); % 

% final x vector
x = zeros(27, 1); % chnaged size from 18 to 27
x(1:6, 1) = [x1; x2];
x(7:27, 1) = xU; 

% final y vector
y = zeros(27, 1);% chnaged size from 15 to 18
y(1:6, 1) = [y1; y2];
y(7:27, 1) = yU;% chnaged size from 1:15 to 1:18

init = [x; y];
end

