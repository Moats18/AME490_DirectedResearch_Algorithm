function [T1, T2, R1, R2] = calc_heli(p, q, theta1, tau1, z, e)
% 
% returns the operators parameterizing a screw isometry generating the
% helical group
%
% inputs:
% (p, q): a pair of integers satisfying the discreteness condition p*tau1 +
% q*tau2 = 0, p*theta1 + q*theta2 = 2*pi.
% theta1: the angle parameterizing the rotation associated with the screw
% isometry used to generate the helical group
% tau1: the translation in R3 associated with the screw isometry used to generate the helical group 
% z: a vector in R3 orthogonal to z parameterizing the origin of the isometry
% e: a unit vector in R3 parameterizing the rotation axis
%
% outputs:
% T1: the translation operator in R3 used in g1
% T2: the translation operator in R3 used in g2
% R1: the rotation operator in SO(3) used in g1
% R2: the rotation operator in SO(3) used in g2

I = eye(3,3); % working in R3

% enforcing the discreteness conditions
tau2 = -p/q*tau1; % p*tau1 = - q*tau2
theta2 = (2*pi - p*theta1)/q; % p*theta1 + q*theta2 = 2*pi

% Rotation matrices where e is the z-axis
R1 = [cos(theta1), -sin(theta1), 0; 
      sin(theta1), cos(theta1), 0;
      0, 0, 1];
R2 = [cos(theta2), -sin(theta2), 0;
      sin(theta2), cos(theta2), 0;
      0, 0, 1];

% Translation operators
T1 = tau1 * e + (I - R1) * z;
T2 = tau2 * e + (I - R2) * z;

end

