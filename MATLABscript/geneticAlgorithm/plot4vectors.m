function plot4vectors(vectors, titles, vis)

% plots four column vectors of 3-dimensional points (x1, y1, z1, x2, y2,
% z2....) separately in 2x2 subplot

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on


for i = 1:length(titles)

subplot(2,2, i)
hold on
plot1vector(vectors(1:27), '', titles{i}, vis);

hold off
end

end
