function plot3dvec(vec, name)

% plots a column vector of 3-dimensional points (x1, y1, z1, x2, y2,
% z2....)

x = zeros(length(vec)/3);
y = x;
z = x;

for i = 1:length(vec)/3
x(i) = vec(3*i-2);
y(i) = vec(3*i-1);
z(i) = vec(3*i);
labels{i} = [name, num2str(i)];
end

plot(x, y);

for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


end

