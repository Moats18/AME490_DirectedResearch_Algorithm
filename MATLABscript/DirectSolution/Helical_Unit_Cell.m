% Initial x-values
p = 6;
q = 4;
x = zeros(1, 3);

tau1 = 0.5;
e1 = [0,0,1]';
z = [0,1,0]';
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

% Visualization setup
figure;
hold on;
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
title('Helical Origami');

% Initialize
numSteps = p*q;
colors = parula(numSteps);
x = x(:); % ensure column vector
panels = zeros(3, numSteps + 1); % store panel origins
panels(:,1) = x;

% Generate transformed panels
for i = 2:numSteps+1
    x_new = R2 * (R1 * panels(:,i-1) + T1) + T2;
    panels(:,i) = x_new;
end


for i = 1:numSteps
    x_curr = panels(:,i);

    g1_x = R1 * x_curr + T1;
    g2_x = R2 * x_curr + T2;
    gg_x = R1 * g2_x + T1;

    % Scatter points
    scatter3(x_curr(1), x_curr(2), x_curr(3), 36, colors(i,:), 'filled');
    scatter3(g1_x(1), g1_x(2), g1_x(3), 36, colors(i,:), 'filled');
    scatter3(g2_x(1), g2_x(2), g2_x(3), 36, colors(i,:), 'filled');
    scatter3(gg_x(1), gg_x(2), gg_x(3), 36, colors(i,:), 'filled');

    % Draw lines
    plot3([x_curr(1) g1_x(1)], [x_curr(2) g1_x(2)], [x_curr(3) g1_x(3)], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot3([x_curr(1) g2_x(1)], [x_curr(2) g2_x(2)], [x_curr(3) g2_x(3)], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot3([g1_x(1) gg_x(1)], [g1_x(2) gg_x(2)], [g1_x(3) gg_x(3)], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot3([g2_x(1) gg_x(1)], [g2_x(2) gg_x(2)], [g2_x(3) gg_x(3)], '-', 'Color', colors(i,:), 'LineWidth', 1.5);
end


