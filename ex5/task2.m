%%

syms a;
A = [0 0 1 0;
     1 2 1 1;
     1 1 0 1;
     1 0 1 2];
C = [0 0 1 0;
     1 1 0 a];

O = [C;C*A];