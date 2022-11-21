clc;
clear all;
syms s f1 f2 f3;

A = [0 1 0; 0 0 1; 1 0 -1];
B = [0; 0; 1];
C = [1 2 0];

Cc = [B A*B A*A*B];
F = [f1 f2 f3];
t = t_beg : 1/f_sampling : t_end - 1/f_sampling;
% simplify(C * inv(s*eye(3) - (A)) * B)
% det(s*eye(3) - A)
expand((s+1)*(s+2)*(2*s+1))