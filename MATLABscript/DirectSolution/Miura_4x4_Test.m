% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Antoine Moats, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 08/18/25
%
% The initial configuration consists of 16 panels in the Miura-Ori
% configuration (4x4 unit cell). This test is interested in the folding of these 
% panels and the resultant energy calculation.

% Test Parameters
test = "shear";
%crease_test = "Miura4x4";
crease_test = "Miura4x4_diagonal2";

% Initial x-values
s = 0.25;

% Idealized ICs
%{
gamma = pi/4;
x1 = [0; 0];
x2 = s*[cos(gamma); sin(gamma)];
x3 = s*[2*cos(gamma); 0];
x4 = x3 + x2;
x5 = 2*x3;
x6 = x5 + s*[0; 1];
x7 = x4 + s*[0; 1];
x8 = s*[2*cos(gamma); 1];
x9 = s*[cos(gamma); sin(gamma)+1];
x10 = s*[0; 1];
x11 = s*[0; 2];
x12 = s*[cos(gamma); sin(gamma)+2];
x13 = s*[2*cos(gamma); 2];
x14 = x13 + x2;
x15 = x6 + x10;
x16 = x15 + x10;
x17 = x14 + x10;
x18 = x13 + x10;
x19 = x12 + x10;
x20 = x11 + x10;
x21 = x20 + x10;
x22 = x19 + x10;
x23 = x18 + x10;
x24 = x17 + x10;
x25 = x16 + x10;
%}

% Randomized ICs
l1R = [1; 0];
l2R = [0; 1];

x1 = [0; 0];
x5 = x1 + l1R;
x21 = x1 + l2R;
x25 = x1 + l1R + l2R;
x3 = (x1+x5)/2;
x11 = (x1+x21)/2;
x15 = (x5+x25)/2;
x13 = (x11+x15)/2;
x23 = (x21+x25)/2;
x2 = (x1+x3)/2;
x4 = (x3+x5)/2;
x6 = (x5+x15)/2;
x8 = (x3+x13)/2;
x7 = (x6+x8)/2;
x10 = (x1+x11)/2;
x9 = (x8+x10)/2;
x12 = (x11+x13)/2;
x14 = (x13+x15)/2;
x16 = (x15+x25)/2;
x18 = (x13+x23)/2;
x17 = (x16+x18)/2;
x20 = (x11+x21)/2;
x19 = (x18+x20)/2;
x22 = (x21+x23)/2;
x24 = (x23+x25)/2;

rng(2, "twister");
r = normrnd(0, 1/32, [50, 1]);
% only perturb eligible DoFs
r(1:2) = 0; % fixed first node
r(41:50) = 0; % top nodes
r(9:12) = 0; % side nodes
r(29:32) = 0; % side nodes

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25];
x = x+r;

% maintain original lattice vectors post-perturbation
x(11) = x(11) + r(19); % shift right side nodes
x(12) = x(12) + r(20);
x(29) = x(29) + r(21);
x(30) = x(30) + r(22);
x(31) = x(31) + r(39);
x(32) = x(32) + r(40);
x(43) = x(43) + r(3); % shift top side nodes
x(44) = x(44) + r(4);
x(45) = x(45) + r(5);
x(46) = x(46) + r(6);
r(47) = x(47) + r(7);
x(48) = x(48) + r(8);

phi = zeros(length(x)/2*3, length(x));
if test == "reference"
    
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
        end
    end
    
elseif test == "axial"
    lambda_1 = 0.50; % axial deformation
    lambda_2 = 1.00;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end

elseif test == "shear"
    %{
    gamma_x = 0;

    % sheared initial x case
    phi_x = zeros(length(x), length(x));
    for i = 1:length(x)
        phi_x(i,i) = 1;
        if mod(i,2) == 1
            phi_x(i,i+1) = gamma_x;
        else
            phi_x(i,i-1) = gamma_x;
        end
    end
    x = phi_x*x;
    %}
    
    gamma = 0.20;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)+1) = gamma;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)-1) = gamma;
        end    
    end

elseif test == "mix"
    gamma = 0.20; % shear deformation
    lambda_1 = 0.50; % axial deformation
    lambda_2 = 0.75;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
            phi(i,i-floor(i/3)+1) = gamma;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
            phi(i,i-floor(i/3)-1) = gamma;
        end    
    end
    % sheared initial x case
    %{
    phi_x = zeros(length(x), length(x));
    for i = 1:length(x)
        phi_x(i,i) = 1;
        if mod(i,2) == 1
            phi_x(i,i+1) = gamma;
        else
            phi_x(i,i-1) = gamma;
        end
    end
    x = phi_x*x;
    %}

end

y = phi*x;

% Bias mountain y's in the z-axis
epsilon = 0.05;
y(18) = y(18)+epsilon;
y(21) = y(21)+epsilon;
y(24) = y(24)+epsilon;
y(27) = y(27)+epsilon;
y(30) = y(30)+epsilon;
y(48) = y(48)+epsilon;
y(51) = y(51)+epsilon;
y(54) = y(54)+epsilon;
y(57) = y(57)+epsilon;
y(60) = y(60)+epsilon;

% x rigidity constraint matrix (11x25) of (2x2)
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2)];
    
% y rigidity constraint matrix (11x25) of (3x3)
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3);
    -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3)];

% populate the vector numbering all of the panels
J = 1:16;

% index set for each panel 
F1 = [1, 2, 9, 10];
F2 = [2, 3, 8, 9];
F3 = [3, 4, 7, 8];
F4 = [4, 5, 6, 7];
F5 = [6, 7, 14, 15];
F6 = [7, 8, 13, 14];
F7 = [8, 9, 12, 13];
F8 = [9, 10, 11, 12];
F9 = [11, 12, 19, 20];
F10 = [12, 13, 18, 19];
F11 = [13, 14, 17, 18];
F12 = [14, 15, 16, 17];
F13 = [16, 17, 24, 25];
F14 = [17, 18, 23, 24];
F15 = [18, 19, 22, 23];
F16 = [19, 20, 21, 22];

if crease_test == "Miura4x4_diagonal1" % middle diagonal crease
    F20 = F16;
    F19 = F15;
    F18 = F14;
    F17 = [17, 24, 25];
    F16 = [16, 17, 25];
    F15 = F12;
    F14 = [13, 14, 17];
    F13 = [13, 17, 18];
    F12 = F10;
    F11 = F9;
    F10 = F8;
    F9 = [9, 12, 13];
    F8 = [8, 9, 13];
    F7 = F6;
    F6 = F5;
    F5 = F4;
    F4 = F3;
    F3 = F2;
    F2 = [1, 2, 9];
    F1 = [1, 9, 10];
    T17 = F17;
    T18 = F18;
    T19 = F19;
    T20 = F20;
    J = 1:20;
end

if crease_test == "Miura4x4_diagonal2" % 2 additional diagonal creases
    F24 = F16;
    F23 = [17; 22; 23];
    F22 = [17; 18; 23];
    F21 = F14;
    F20 = [17; 24; 25];
    F19 = [16; 17; 25];
    F18 = F12;
    F17 = [13; 14; 17];
    F16 = [13; 17; 18];
    F15 = F10;
    F14 = [11; 12; 19];
    F13 = [11; 19; 20];
    F12 = F8;
    F11 = [9; 12; 13];
    F10 = [8; 9; 13];
    F9 = F6;
    F8 = [7; 14; 15];
    F7 = [6; 7; 15];
    F6 = F4;
    F5 = [3; 4; 7];
    F4 = [3; 7; 8];
    F3 = F2;
    F2 = [1; 2; 9];
    F1 = [1; 9; 10];
    T17 = F17;
    T18 = F18;
    T19 = F19;
    T20 = F20;
    T21 = F21;
    T22 = F22;
    T23 = F23;
    T24 = F24;
    J = 1:24;
end

T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;
T5 = F5;
T6 = F6;
T7 = F7;
T8 = F8;
T9 = F9;
T10 = F10;
T11 = F11;
T12 = F12;
T13 = F13;
T14 = F14;
T15 = F15;
T16 = F16;

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end
%{
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Tj(1, :, 3) = T3;
Tj(1, :, 4) = T4;
Tj(1, :, 5) = T5;
Tj(1, :, 6) = T6;
Tj(1, :, 7) = T7;
Tj(1, :, 8) = T8;
Tj(1, :, 9) = T9;
Tj(1, :, 10) = T10;
Tj(1, :, 11) = T11;
Tj(1, :, 12) = T12;
Tj(1, :, 13) = T13;
Tj(1, :, 14) = T14;
Tj(1, :, 15) = T15;
Tj(1, :, 16) = T16;
%}
Fj = Tj;

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, crease_test);
