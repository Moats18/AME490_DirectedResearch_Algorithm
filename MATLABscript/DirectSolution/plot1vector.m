function plot1vector(vec, name, title1, vis, test)

% plots a column vector of 3-dimensional points (x, y, z1, x2, y2,
% z2....)

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on

gap = 0.4;
vec = cell2mat(vec);
x = zeros(length(vec)/2, 1);
y = zeros(length(vec)/2, 1);

for i = 1:length(vec)/2
x(i) = vec(2*i-1);
y(i) = vec(2*i);
labels{i} = [name, num2str(i)];
end

figure
title(title1)
xlabel('x');
ylabel('y');

hold on
xlim([min(x) - gap, max(x) + gap]);
ylim([min(y) - gap, max(y) + gap]);
plot(x, y);

for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


% Use the plot function to draw the horizontal line
if test == "Miura"
    plot([x(1), x(6)], [y(1), y(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(5)], [y(2), y(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(5), x(8)], [y(5), y(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(4), x(9)], [y(4), y(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Folding"
    plot([x(1), x(6)], [y(1), y(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(5)], [y(2), y(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Rotating Squares"
    plot([x(1), x(2)], [y(1), y(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(12)], [y(1), y(12)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(12), x(3)], [y(12), y(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(12), x(9)], [y(12), y(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(6)], [y(3), y(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(7), x(8)], [y(7), y(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(6), x(9)], [y(6), y(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

%visualization of the lattice vectors 
%{
if vis == 1
plot([x(1) - 0.5*gap, x(1) - 0.5*gap], [y(10), y(2)], 'r-', 'LineWidth', 0.5);
plot([x(1) - 0.3, x(1) - 0.1], [y(10), y(10)], 'r-', 'LineWidth', 0.5);
plot([x(1) - 0.3, x(1) - 0.1], [y(2), y(2)], 'r-', 'LineWidth', 0.5);
h = text(x(1) -0.04, (y(11)+y(1))/1.5, "E1 = " + num2str(y(10) - y(2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
set(h,'Rotation',90);

plot([x(1), x(5)], [y(1)- 0.5*gap, y(1)- 0.5*gap], 'r-', 'LineWidth', 0.5); 
plot([x(1), x(1)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
plot([x(5), x(5)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
text((x(1)+x(5))/1.5, y(1)-0.5*gap, "E2 = " + num2str(x(5) - x(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
end
%}
hold off

end

