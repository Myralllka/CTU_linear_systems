syms s

A = [1 4 0; 0 -1 0; 0 2 -3];
B = [0; -1; -1];
C = [-1 -1 0];
D = [0];
I = [1 0 0; 0 1 0; 0 0 1];

Ps = [s*I - A, B; -C, D]
smithForm(Ps)
rank(Ps)

H = C * (inv(s*I - A)) *B;

S = MNsmithForm([s*I - A;C]);

[z,p,k] = ss2zp(A,B,C,D);

Ct = [C; C*A; C*A*A];
null(Ct);
Q = [0 1 0 ; 
     1 0 0 ; 
     0 0 1];

Ac = inv(Q)*A*Q;

% -----------------------------

Nr = [s-1, 2; s+1 0];
Dr = [s+1, 0; 0, s^2-1];

rank([Nr, Dr]);
MNsmithForm([Nr; Dr]);

Gr = [-2, 2; (s+1), 0];
(Nr * inv(Gr)) * inv(Dr * inv(Gr));
MNsmithForm([Nr*inv(Gs); Dr*inv(Gs)]);
simplify(Nr * inv(Gr));
simplify(Dr * inv(Gr));

function [SA, invFact, D] = MNsmithForm(A)
    % from https://www.mathworks.com/matlabcentral/fileexchange/60410-mnsmithform-a
    row = size(A,1);
    col = size(A,2);
    n = min(row,col);
    minors = cell(1,n);
    D0 = 1;
    D0 = sym(D0);
    D = sym(NaN(1,n));
    invFact = sym(NaN(1,n));
    for i = 1:n;
        rowindex = false(1,row);
        rowindex(1:i) = true;
        rowperms = unique(perms(rowindex),'rows');
        colindex = false(1,col);
        colindex(1:i) = true;
        colperms = unique(perms(colindex),'rows');
        rownum = size(rowperms,1);
        colnum = size(colperms,1);
        minors{i} = sym(NaN(rownum,colnum));
        for j=1:rownum;
            for k=1:colnum;
                Atmp = A;
                Atmp = Atmp(rowperms(j,:),:);
                Atmp = Atmp(:,colperms(k,:));
                minors{i}(j,k) = det(Atmp);
            end
        end
        rowlen = rownum*colnum; %(row - (i-1))*(col - (i-1));
        minors{i} = reshape(minors{i},1,rowlen);
        minors{i}(minors{i} == 0) = [];
        D(i) = gcd(minors{i});
        if i == 1
            invFact(i) = D(i)/D0;
        else
            invFact(i) = D(i)/D(i-1);
        end   
    end
    SA = diag(invFact);
    if row>col
        zerorows = zeros(row-col,col);
        SA = [SA;zerorows];
    elseif col>row
        zerocols = zeros(row,col-row);
        SA = [SA, zerocols];
    end
end