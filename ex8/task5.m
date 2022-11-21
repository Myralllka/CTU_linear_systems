clc;
clear all;

syms om s p1 p2 p3 p4 r;

R = [r];
assume(r>0);
assume(om>0);
A = [0 1; pi/2 0];
B = [0; 1];
C = [1 1];
P = [p1 p2; p3 p4];
Q = [1 1; 1 1];

res = A.' * P + P*A - P*B*inv(R)*B.'*P + Q;
[rp1, rp2, rp3, rp4] = solve([res(1) == 0, res(2) == 0, res(3) == 0, res(4) == 0, p2==p3, det(res(1)) >= 0, det(P) >= 0], [p1, p2, p3, p4]);
simplify(rp1);
simplify(rp2);
simplify(rp3);
simplify(rp4);
P1 = [rp1(1), rp2(1); rp3(1), rp4(1)];
% P2 = [rp1(2), rp2(2); rp3(2), rp4(2)];
% P3 = [rp1(3), rp2(3); rp3(3), rp4(3)];
% P4 = [rp1(4), rp2(4); rp3(4), rp4(4)];
F = simplify(-(inv(R) + B.' * P1));

% sys = ss(A + B*F, B, C, [0])
X = [];
Y = [];

figure(1)
hold on
counter = 1;
for t = 100:-0.1:1
    e = eig(subs(A+B*F, r, t));
    
    X(counter) = real(e(1));
    Y(counter) = imag(e(1));
    counter = counter + 1;
    X(counter) = real(e(2));
    Y(counter) = imag(e(2));
    counter = counter + 1;
end

scatter3(X, Y, [counter:-0.1:1])

xlabel('Re');
ylabel('Im');
zlabel("r");

% [X,K,L] = icare(A,B,Q,R,[],eye(2),eye(2));
