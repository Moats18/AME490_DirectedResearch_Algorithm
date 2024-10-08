function subplot3dvec(vec1, name1, vec2, name2)

% plots two column vectors of 3-dimensional points (x1, y1, z1, x2, y2,
% z2....) separely in a horizontal subplot

x1 = zeros(length(vec1)/3, 1);
y1 = x1;
z1 = x1;

x2 = x1;
y2 = y1;
z2 = z1;

for i = 1:length(vec1)/3
x1(i) = vec1(3*i-2);
y1(i) = vec1(3*i-1);
z1(i) = vec1(3*i);
labels1{i} = [name1, num2str(i)];
end

for i = 1:length(vec1)/3
x2(i) = vec2(3*i-2);
y2(i) = vec2(3*i-1);
z2(i) = vec2(3*i);
labels2{i} = [name2, num2str(i)];
end

subplot(1, 2, 1);
hold on
title('Initial Configuration')
xlabel('x');
ylabel('y');
plot(x1, y1);
for i = 1:length(labels1)
    text(x1(i), y1(i), labels1{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


% Use the plot function to draw the horizontal line
plot([x1(1), x1(6)], [y1(1), y1(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
plot([x1(2), x1(5)], [y1(2), y1(5)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

hold off


subplot(1, 2, 2);
hold on
title('Final Configuration')
xlabel('x');
ylabel('y');
plot(x2, y2);
for i = 1:length(labels2)
    text(x2(i), y2(i), labels2{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


plot([x2(1), x2(6)], [y2(1), y2(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
plot([x2(2), x2(5)], [y2(2), y2(5)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

hold off


end

