A = [1 4 0; 2 -1 0; 0 0 1];
B = [0 0; 1 0; -1 1];
C = [B, A*B, A*A*B];
% rank(C)
% eig(A)
rank(C(:, 1:3));
Cn = [0 4 0; 1 -1 0; -1 -1 1];
inv(Cn);
P = [1/4 0 0; 1/4 1 0; 1/2 1 1];
Pi = inv(P);

Ac = P*A*Pi;
Bc = P*B;

syms s;
I = eye(3);

Fc1 = [-9 0 1; -9 -3 -4];

K = [ 0  1 0; 
      -1 -2 0;
      0 0 -1];

Fc2 = K(2:3, :) - [9 0 0; 8 0 1];

F1 = Fc1*P
F2 = Fc2*P

det(s*I - (A + B*F1))
jordan((A + B*F1))
det(s*I - (A + B*F2))
jordan((A + B*F2))


% det(s*eye(3) - K)
% jordan(K)