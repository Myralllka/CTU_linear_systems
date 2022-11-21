A = [2 1 0 0; 0 2 0 0; 0 0 -1 0; 0 0 0 -1];
B = [0; 1; 1; 1];
C = [B, A*B, A*A*B, A*A*A*B];
% eig(A)
rank(C);
Q = [0 1 12 0; 1 2 8 0; 1 -1 -1 0; 1 -1 -1 1];
% rank(Q)
Ac = inv(Q) * A * Q;
Bc = inv(Q) * B;
A1 = [0 4/3 0; 1 0 -4; 0 1/3 3];
eig(A1);
syms s f1 f2 f3 f4;
F = [f1 f2 f3 f4];
I = eye(4);
simplify(det(s*I - (Ac + Bc*F)))
