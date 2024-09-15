function plot4vectors3D(vectors, titles, vis)

% plots four column vectors of 3-dimensional points (x1, y1, z1, x2, y2,
% z2....) separately in 2x2 subplot

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on


for i = 1:length(titles)

plot1vector3D(vectors(:, i), '', titles{i}, vis);

end

end
