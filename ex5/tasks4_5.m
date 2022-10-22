%%

syms om  s;

A = [0 1 0 0;
    3*om^2 0 0 2*om;
    0 0 0 1;
    0 -2*om 0 0];

B = [0; 0; 0; 1];

Q = [0 1 0 0;
     1 0 -om^2 0;
     0 0 -2*om 0;
     0 -2*om 0 1];
I = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

C = [B A*B A*A*B A*A*A*B];
tmp = inv(C);
q = tmp(4, :);
P = [q;q*A;q*A*A;q*A*A*A];
Ac = P*A*inv(P);

det(s*I - A);
A = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     0 0 -om^2 0];

C = [B A*B A*A*B A*A*A*B];

rank(C)