syms s

A = [1 0; 0 -1];
B = [0 1; 2 1];
C = [1 -1; 0 0];
D = [1 0; 1 0];
I = [1 0; 0 1];

H = C * (inv(s*I - A)) *B + D

% [z,p,k] = ss2zp(A,B,C,D);
% eig(A);