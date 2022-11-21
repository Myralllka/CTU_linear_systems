syms s om;
A = [0 1 0 0; 9*om^2 0 0 2*om; 0 0 0 1; 0 -2*om -4*om^2 0];
B = [0; 0; 0; -1];

simplify(expand((s+3*om)*(s+4*om)*(s+3*om+3*i*om)*(s+3*om-3*i*om)));
Ci = inv([B, A*B, A*A*B, A*A*A*B]);
q = Ci(end, :);
P = [q;q*A; q*A*A; q*A*A*A]
Ac = P*A*inv(P);
Bc = inv(P) * B;
syms f1 f2 f3 f4;
F = [f1 f2 f3 f4];
% det(s*eye(4) - (Ac + B*F))
a1 = s + 3*om;
a2 = s + 4*om;
a3 = s + (1+i)*3*om;
a4 = s + (1-i)*3*om;

des = expand(a1*a2*a3*a4);

Fc = [-252*om^4 -198*om^3 -73*om^2 -13*om];
F = Fc*P;

eig(A+B*F)
eig(Ac+Bc*Fc)
% roots([1 13*om 72*om^2 198*om^3 216*om^4])