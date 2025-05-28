% Initial x-values
p = 6;
q = 4;
x = zeros(p*q, 3);

index = 1;
for i = 1:p
    for j = 1:q
        x(index, 1) = i - 1;
        x(index, 2) = j - 1;
        x(index, 3) = 0;
        index = index + 1;
    end
end

tau1 = 0.5;
e1 = [0,0,1]';
z = [1,0,0]';
theta1 = pi/4;

% discreteness condition 
tau2 = -p/q*tau1;
theta2 = 2*pi/q - p/q*theta1;
I = eye(3,3);

% Rotation matrices
R1 = [cos(theta1), -sin(theta1), 0; sin(theta1), cos(theta1), 0; 0, 0, 1];
R2 = [cos(theta2), -sin(theta2), 0; sin(theta2), cos(theta2), 0; 0, 0, 1];

% Translation vectors
T1 = tau1 * e1 + (I - R1) * z;
T2 = tau2 * e1 + (I - R2) * z;

% Transformation
n = length(x);
x_current = x;
numSteps = 20;

% Visualization setup
figure;
hold on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
colors = parula(numSteps);

for i = 1:numSteps
    % Plot the current step
    scatter3(x_current(1,1), x_current(1,2), x_current(1,3), 36, colors(i,:), 'filled');

    scatter3(x_current(6,1), x_current(6,2), x_current(6,3), 36, colors(i,:), 'filled');

    scatter3(x_current(2,1), x_current(2,2), x_current(2,3), 36, colors(i,:), 'filled');
    scatter3(x_current(5,1), x_current(5,2), x_current(5,3), 36, colors(i,:), 'filled');

   % Store previous state
    x_prev = x_current;

    x_temp = R1 * x_current' + T1;    
    x_next = R2 * x_temp + T2;         
    x_next = x_next';                 

     % Draw lines between panel points

        plot3([x_current(1,1), x_current(2,1)], ...
              [x_current(1,2), x_current(2,2)], ...
              [x_current(1,3), x_current(2,3)], ...
              'Color', colors(i,:), 'LineWidth', 1);

        plot3([x_current(5,1), x_current(6,1)], ...
              [x_current(5,2), x_current(6,2)], ...
              [x_current(5,3), x_current(6,3)], ...
              'Color', colors(i,:), 'LineWidth', 1);

      plot3([x_current(2,1), x_current(6,1)], ...
              [x_current(2,2), x_current(6,2)], ...
              [x_current(2,3), x_current(6,3)], ...
              'Color', colors(i,:), 'LineWidth', 1);

       plot3([x_current(1,1), x_current(5,1)], ...
              [x_current(1,2), x_current(5,2)], ...
              [x_current(1,3), x_current(5,3)], ...
              'Color', colors(i,:), 'LineWidth', 1);
    
    % Draw lines between corresponding layers 
   
    for j = [1, 2, 5, 6]
        plot3([x_prev(j,1), x_next(j,1)], ...
              [x_prev(j,2), x_next(j,2)], ...
              [x_prev(j,3), x_next(j,3)], ...
              'Color', colors(i,:), 'LineWidth', 1);
    end
    %}
    x_current = x_next;
end





