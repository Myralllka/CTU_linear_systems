syms s;
P1 = [s^2 + s, -s;
      -s^2-1, s^2];
P2 = [s, 0;
      -s-1, 1];

% U = gcd(P1, P2)
inv(P2) * P1