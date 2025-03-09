function [U] = Unitary(p)
%[U] = UNITARY(M)
%Calculates a random complex unitary matrix, size p.modes*p.modes

if isnumeric(p)  && p <= 0
    U =  'complex unitary';
    return;
end
U = (randn(p.modes) + 1i*randn(p.modes))/(sqrt(2));
[U,R] = qr(U);
D = diag(R);
ph = D./abs(D);
U = U*diag(ph);
end