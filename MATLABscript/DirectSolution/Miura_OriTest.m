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
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.

test = "mix";

% Initial x-values
s = 0.5;

gamma = pi/4;
x1 = [0; 0];
x2 = s*[cos(gamma); sin(gamma)];
x3 = s*[2*cos(gamma); 0];
x4 = s*[2*cos(gamma); 1];
x5 = s*[cos(gamma); sin(gamma)+1];
x6 = s*[0; 1];
x7 = s*[0; 2];
x8 = s*[cos(gamma); sin(gamma)+2];
x9 = s*[2*cos(gamma); 2];

%{
gamma = pi/6;
beta = pi*7/18;
x1 = [0; 0];
x2 = s*[cos(gamma); sin(gamma)];
x3 = s*[2*cos(gamma); 0];
x6 = s*[cos(beta); sin(beta)];
x7 = 2*s*[cos(beta); sin(beta)];
x9 = x3+x7;
x4 = (x3+x9)/2;
x8 = x2 + x9-x3;
x5 = (x2+x8)/2;
%}
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

l1R = x3 - x1;
l2R = x7 - x1;

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
    lambda_1 = 1.00; % axial deformation
    lambda_2 = 1.00;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2ff
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end

elseif test == "shear"
    gamma = 0.10;
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
    lambda_1 = 1; % axial deformation
    lambda_2 = sqrt(2)/2;
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

% Bias y's (y4 y5 y6) in the z-axis
y(12) = y(12)+0.1;
y(15) = y(15)+0.1;
y(18) = y(18)+0.1;

% known Miura solution
%{
%gamma = 0;
y1 = [0; 0; 0];
y2 = s*[cos(gamma); sin(gamma); 0];
y3 = s*[2*cos(gamma); 0; 0];
y4 = s*[2*cos(gamma); 1; 0] + [0; 0; 1];
y5 = s*[cos(gamma); sin(gamma)+1; 0] + [0; 0; 1];
y6 = s*[0; 1; 0] + [0; 0; 1];
y7 = s*[0; 2; 0];
y8 = s*[cos(gamma); sin(gamma)+2; 0];
y9 = s*[2*cos(gamma); 2; 0];
y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];

% First, some useful parameters of the Miura-Ori fold
gamma = pi/2 - gamma;
theta = pi/6;
lambda = 1;
H = s*sin(theta)*sin(gamma);
S = s*cos(theta)*tan(gamma)/sqrt(1+cos(theta)^2*tan(gamma)^2);
L = s*sqrt(1-sin(theta)^2*sin(gamma)^2);
V = s*1/sqrt(1+cos(theta)^2*tan(gamma)^2);
eta = atan(cos(theta)*tan(gamma));
psi = asin(sin(theta)*sin(gamma));
phi = asin(sin(eta)/sin(gamma));
o = H/tan(theta);

% Initial y-values
y1 = lambda*[0; 0; 0];
y2 = lambda*[S; s*cos(eta); 0];
y3 = lambda*[2*S; 0; 0];
y4 = lambda*[2*S; L; H];
y5 = lambda*[S; 2*L; H];
y6 = lambda*[0; L; H];
y7 = lambda*[0; 2*L; 0];
y8 = lambda*[S; 2*L+V; 0];
y9 = lambda*[2*S; 2*L; 0];
y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];
%}


% x rigidity constraint matrix (7x9) of (2x2)
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);];
    
% y rigidity constraint matrix (7x9) of (3x3)
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);];

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
J = [1, 2, 3, 4];

% index set for each panel 
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
F3 = [5, 6, 7, 8];
F4 = [4, 5, 8, 9];
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;

% 3D array containing index set of x coordinates for panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Tj(1, :, 3) = T3;
Tj(1, :, 4) = T4;
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

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");
