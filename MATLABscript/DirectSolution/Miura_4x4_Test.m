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
% Updated Date: 07/22/25.
%
% The initial configuration consists of 16 panels in the Miura-Ori
% configuration (4x4 unit cell). This test is interested in the folding of these 
% panels and the resultant energy calculation.

test = "shear";

% Initial x-values
s = 0.25;

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

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25];

l1R = x5 - x1;
l2R = x21 - x1;

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
    gamma = 0.05;

    % sheared initial x case
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
    
    gamma = 0.05;
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
epsilon = 0.5;
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

% center point buckle case
%{
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);];
    
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);];
%}

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

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura4x4");
