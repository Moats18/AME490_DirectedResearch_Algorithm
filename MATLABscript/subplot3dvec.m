function subplot3dvec(vec1, name1, vec2, name2)

x1 = zeros(length(vec1)/3);
y1 = x1;
z1 = x1;

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
hold off


end

