function [T1, T2, R1, R2] = calc_heli(q, p, theta1, tau1, z, e1)

tau2 = -p/q*tau1;
theta2 = 2*pi/q - p/q*theta1;
I = eye(3,3);

% Rotation matrices where e1 is the z-axis
R1 = [cos(theta1), -sin(theta1), 0; sin(theta1), cos(theta1), 0; 0, 0, 1];
R2 = [cos(theta2), -sin(theta2), 0; sin(theta2), cos(theta2), 0; 0, 0, 1];

% Translation vectors
T1 = tau1 * e1 + (I - R1) * z;
T2 = tau2 * e1 + (I - R2) * z;


end

