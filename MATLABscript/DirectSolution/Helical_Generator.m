function [points, R1, T1, R2, T2] = Helical_Generator(p, q, tau1, e1, z, theta1)

x = zeros(1, 3);

[T1, T2, R1, R2] = calc_heli(q, p, theta1, tau1, z, e1);

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
numSteps = p*(q-1)/2;
colors = parula(numSteps);
x = x(:); % ensure column vector
panels = zeros(3, numSteps + 1); % store panel origins
panels(:,1) = x;
points_all = [];

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

    % Store all generated points
    points_all = [points_all; x_curr'; g1_x'; g2_x'; gg_x'];

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

rounded_points = round(points_all, 1);
points = unique(rounded_points, 'rows');

end

