clc;
clear all;

syms s a;
A = [0 0; 1 -1];
B = [0 1; -1 1];
C = [0 1];
D = [1 0];
I = [1 0; 0 1];
H = C*inv(s*I - A)*B + D;
simplify(H)

A = [0 0 0; 
     1 0 0;
     0 1 -1];
B = [0, 0; 
     B];
C = [0 0 1];
D = [1 0];
I = eye(3);
Oo = [C;C*A;C*A*A];
Cc = [B, A*B, A*A*B];
H = C*inv(s*I - A)*B + D;
simplify(H)

rank(Oo)
rank(Cc)
% H = [s+1 a; s+2 s+2];
% hh = smithForm(H)
% simplify(hh/(s*(s+2)));

% N = [(s+1), a; s+2, s+2];
% R = [s*(s+2), 0; 0, s*(s+2)];

% simplify(N*inv(R));