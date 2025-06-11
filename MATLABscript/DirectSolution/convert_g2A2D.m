function A = convert_g2A2D(R1, R2, T1, T2)

% determining the constraint matrix A from the g = R(x) + C 
% so that Ay = b using a 2x2 unit cell (nine points) 
% added x6 -> x4 constraint (incompatible with Kirigami and Rotating
% Squares)
% 
% Labeling System
% 7 -> 8 -> 9 
% 6 <- 5 <- 4 
% 1 -> 2 -> 3

% adding the constraint of the first node being fixed
A = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     R1,    zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), R1, zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), R1, zeros(2), -eye(2);
     R2,    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2);
     zeros(2), R2,    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), zeros(2), R2,    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2)]; % x3 to x9 is a redundant constraint

end


