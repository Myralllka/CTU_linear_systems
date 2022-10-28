syms s

A = [2 0; 1 1];
B = [1; 0];
C = [1 1];
D = [0];
I = [1 0; 0 1];

H = C * (inv(s*I - A)) *B + D
simplify(H - (s/((s-1)*(s-2))))
[z,p,k] = ss2zp(A,B,C,D);
eig(A);