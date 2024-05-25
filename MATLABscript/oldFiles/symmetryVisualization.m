function symmetryVisualization(vec1, name)

x1 = zeros(length(vec1)/3, 1);
y1 = x1;
z1 = x1;

for i = 1:length(vec1)/3
x1(i) = vec1(3*i-2);
y1(i) = vec1(3*i-1);
z1(i) = vec1(3*i);
labels1{i} = [name, num2str(i)];
end

figure
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

plot([x1(1), x1(1)], [y1(1), y1(3)], 'r-', 'LineWidth', 0.5); 
plot([x1(1), x1(6)], [y1(1), y1(1)], 'r-', 'LineWidth', 0.5);

hold off


end

