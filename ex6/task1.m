
syms s;
A = [(s-1)*(s+1) 2*(s-2); 
     s*(s-2) 0];

[U, V, S] = smithForm(A)
minor = @(i,j,A)det(A(setdiff([1:end],[i]),setdiff([1:end],[j])));

% minor(2, 2, A)
det(A)