function plot1vector3D(vec, name, title1, vis, test)

% plots a column vector of 3-dimensional points (x, y, z1, x2, y2,
% z2....)

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on

gap = 0.4;
vec = cell2mat(vec);
x = zeros(length(vec)/3, 1);
y = x;
z = x;

for i = 1:length(vec)/3
x(i) = vec(3*i-2);
y(i) = vec(3*i-1);
z(i) = vec(3*i);
labels{i} = [name, num2str(i)];
end


figure
xlim([min(x) - gap, max(x) + gap]);
ylim([min(y) - gap, max(y) + gap]);
plot3(x, y, z);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

hold on

for i = 1:length(labels)
    text(x(i), y(i), z(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

title(title1)

%Use the plot function to draw lines

if test == "Miura"
    plot3([x(1), x(6)], [y(1), y(6)], [z(1), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(5)], [y(2), y(5)], [z(2), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(5), x(8)], [y(5), y(8)], [z(5), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(4), x(9)], [y(4), y(9)], [z(4), z(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Folding"
    plot3([x(1), x(6)], [y(1), y(6)], [z(1), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(5)], [y(2), y(5)], [z(2), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Rotating Squares"
    plot3([x(1), x(2)], [y(1), y(2)], [z(1), z(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(1)], [y(12), y(1)], [z(12), z(1)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(3)], [y(12), y(3)], [z(12), z(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(9)], [y(12), y(9)], [z(12), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3), x(6)], [y(3), y(6)], [z(3), z(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(7), x(8)], [y(7), y(8)], [z(7), z(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(6), x(9)], [y(6), y(9)], [z(6), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "1DBar"
    plot3([x(1), x(2)], [y(1), y(2)], [z(1), z(2)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(3)], [y(2), y(3)], [z(2), z(3)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

%{
for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

%visualization of the lattice vectors 
if vis == 1
plot([x(1) - 0.5*gap, x(1) - 0.5*gap], [y(1), y(3)], 'r-', 'LineWidth', 0.5);
plot([x(1) - 0.3, x(1) - 0.1], [y(3), y(3)], 'r-', 'LineWidth', 0.5);
plot([x(1) - 0.3, x(1) - 0.1], [y(1), y(1)], 'r-', 'LineWidth', 0.5);
h = text(x(1) -0.04, (y(3)+y(1))/1.5, "E1 = " + num2str(y(3) - y(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
set(h,'Rotation',90);

plot([x(1), x(6)], [y(1)- 0.5*gap, y(1)- 0.5*gap], 'r-', 'LineWidth', 0.5); 
plot([x(1), x(1)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
plot([x(6), x(6)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
text((x(1)+x(6))/1.5, y(1)-0.5*gap, "E2 = " + num2str(x(6) - x(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
end
%}
hold off
end
