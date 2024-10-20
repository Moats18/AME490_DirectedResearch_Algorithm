
x1 = [0; 0; 0];
x2 = [1; 0; 0];
x3 = [1; 1; 0];
x4 = [0; 1; 0];
x = [x1; x2; x3; x4];

theta = pi/6;
R = [1  0           0;
     0  cos(theta)  -sin(theta);
     0  sin(theta)  cos(theta)];

disp("Given rotation matrix: ")
disp(R)

y1 = R*x1;
y2 = R*x2;
y3 = R*x3;
y4 = R*x4;
y = [y1; y2; y3; y4];

[c, r] = centerOfPanel([1, 2, 3, 4], x);

for i = 1:4
       V = r(3*i-2:3*i) + c - y(3*i-2:3*i, 1);
       T = c - r(3*i-2:3*i) - y(3*i-2:3*i);
       Bi{i} =  [0       V(1)    V(2)    V(3);
                T(1)    0       -V(3)   V(2);
                T(2)    V(3)    0       -V(1);
                T(3)    -V(2)   V(1)    0];
end

B = zeros(4);

    for i = 1:4
    B = B + Bi{i}'*Bi{i};
    end

    % rounding the small numeric values to zero
    roundedB = B;
    roundedB(abs(roundedB)<1e-3)=0;

    [V ,D] = eig(roundedB);
    roundedD = D;
    roundedD(abs(roundedD)<1e-3)=0;
    maxEig = max(roundedD, [], "all");
    [~, col] = find(roundedD == maxEig);

    if length(col) > 1
        eigenVal = col(1);
        disp("Multiple max eigenvalues!")
    else
        eigenVal = col;
        disp("One max eigenvalue!")
    end

    pj = V(:, eigenVal);
    v = pj(2:4);
    Rnew = quat2rotm(pj');

for i = 1:4
pj = V(:, i);
Rnew = quat2rotm(pj');

e = 0;

if i ==  eigenVal
disp("Optimized rotation matrix: ")
disp(Rnew);
else 
fprintf('Rotation matrix of %d quaternion\n', round(i));
disp(Rnew);

end

for j = 1:4
    e = e + norm(y(3*j-2:3*j, 1) - Rnew*r(3*j-2:3*j, 1) - c)^2;
end

disp(e);

end

Fj = [1,2,3,4];
Tj = Fj;
J = 1;


%Ropt = rotationMin_new(x, y, Fj, Tj, J);
%disp(Ropt{1})

%{
for i = 1:length(x)/3
xx(i) = x(3*i-2);
xy(i) = x(3*i-1);
xz(i) = x(3*i);
end

for i = 1:length(y)/3
yx(i) = y(3*i-2);
yy(i) = y(3*i-1);
yz(i) = y(3*i);
end

figure;
hold on;
plot3(xx, xy, xz);
plot3(yx, yy, yz);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

%}

