clc;
clear all;

syms s m1 m2 m3 m4 b1 b5 k1 k2 k3 k4 k5;

M = [m1 0 0 0; 
     0 m2 0 0; 
     0 0 m3 0; 
     0 0 0 m4];
B = [b1 0 0 0; 
     0 0 0 0; 
     0 0 0 0; 
     0 0 0 b5];
K = [k1+k2 -k2 0 0; 
    -k2 k2+k3 -k3 0; 
    0 -k3 k3+k4 -k4; 
    0 0 -k4 k4+k5];
G = [0; 1; 0; 0];
P = [1 0 0 0];

M = subs(M, m1, 1);
M = subs(M, m4, 1);
M = subs(M, m2, 2);
M = subs(M, m3, 2);
K = subs(K, k1, 1);
K = subs(K, k4, 1);
K = subs(K, k2, 2);
K = subs(K, k3, 2);
K = subs(K, k5, 4);
B = subs(B, b1, 1);
B = subs(B, b5, 2);

I = eye(4);
O = zeros(4, 4);

Ass = eval([O I; -inv(M)*K -inv(M)*B]);
Bss = eval([0; 0; 0; 0; inv(M) * G]);
Css = [P 0 0 0 0];
Dss = [0];

% ----------------------------------------

t = 0:0.04:100;  % 201 points
u = max(0,min(t-1,1));

figure(1);
%
subplot(2, 2, 1);
sys = ss(Ass, Bss, Css, [0]);
pzmap(sys)
title('original system zeros poles');
grid on
%
subplot(2, 2, 3);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 8);
sys = ss(a, b, c, [0]);
pzmap(sys)
title('original system balanced realization zeros and poles, order 8');
grid on
%
subplot(2, 2, 2);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 6);
sys = ss(a, b, c, [0]);
pzmap(sys)
title('original system balanced realization zeros and poles, order 6');
grid on
% 
subplot(2, 2, 4);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 4);
sys = ss(a, b, c, [0]);
pzmap(sys)
title('original system balanced realization zeros and poles, order 4');
grid on

% ==============================================================================
figure(2);
%
subplot(2, 2, 1);
sys = ss(Ass, Bss, Css, [0]);
lsim(sys,u,t);
title('original system simulation');
grid on
%
subplot(2, 2, 3);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 8);
sys = ss(a, b, c, [0]);
lsim(sys,u,t);
title('original system balanced realization simulation, order 8');
grid on
%
subplot(2, 2, 2);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 6);
sys = ss(a, b, c, [0]);
lsim(sys,u,t)
title('original system balanced realization simulation, order 6');
grid on
% 
subplot(2, 2, 4);
[a, b, c] = balanced_realisation_particular(Ass, Bss, Css, Dss, 4);
sys = ss(a, b, c, [0]);
lsim(sys,u,t)
title('original system balanced realization simulation, order 4');
grid on

% ==============================================================================
figure(3);
H = Css * inv(s*eye(8) - Ass) * Bss + Dss;
% [a, b, c] = balanced_realisation_neg(H, 4);
[h] = balanced_realisation_neg(H, 4);
h = simplifyFraction(simplify(h))
[n, d] = numden(h);
n = coeffs(n);
d = coeffs(d);
d = fliplr(d);
sys = tf(eval(n), eval(d));

% sys = ss(a, b, c, [0]);
lsim(sys,u,t)
% title('order 4 without fastest poles');
grid on


function [A, B, C] = balanced_realisation_particular(Ass, Bss, Css, Dss, m)
    syms s;
    H = Css * inv(s*eye(8) - Ass) * Bss + Dss;
    [~, r] = numden(H);
    r = polynomialDegree(r);
    
    A_pows = zeros(size(Ass, 1), size(Ass, 2), r*2);
    A_pows(:, :, 1) = eye(8);
    for i=2:r*2
       A_pows(:, :, i) = A_pows(:, :, i-1) * Ass;
    end
    T = zeros(r, r);
    That = zeros(r, r);
    
    for i=1:r
        for j=1:r
            That(i, j) = Css * A_pows(:, :, i+j) * Bss;
            T(i, j) = Css * A_pows(:, :, i+j-1) * Bss;
        end
    end

    [K, S, Lt] = svd(T);
    L = Lt.';
    K = K(:, 1:m);
    L = L(1:m, :);
    S = S(1:m, 1:m);
    Vplus_8 = S^(-1/2) * K.';
    Uplus_8 = L.' * S^(-1/2);
    V_8 = K * S^(1/2);
    U_8 = S^(1/2) * L;

    A = Vplus_8 * That * Uplus_8;
    B = U_8(:, 1); % because H 1x1
    C = V_8(1, :); % because H 1x1
    
end

% function [A, B, C] = balanced_realisation_neg(H, m)
function [H_new] = balanced_realisation_neg(H, m)

    % Ok/ I have spent 14 hours on it and I have no idea how to do it correctly.

    syms s;

    a = sort(poles(H), 'descend');
    
    dc_original = limit(H, s, 0); % DC gain original 

    H = 1/(s-a(1));
    
    % ignore 4 zeros and 4 fastest poles
    for i=2:m
        H = H / (s-a(i));
    end
    dc_new = limit(H, s, 0); % DC gain new
    H = H * dc_original / dc_new; % to preserve DC gain
    
    H_new = H;
   
end