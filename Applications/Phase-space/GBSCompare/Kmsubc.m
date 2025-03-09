function C = Kmsubc(s)
% C = KMSUBC(a,p); Generates all a subset of all possible combinations of
% output mode click correlation moments for sequential modes. This is
% done for graphical simplificity. Assumes output unchanged from input. 
% eg, identity matrix or uniform thermal+unitary

n = (sinh(s.r').*s.t').^2;
m = ((1-s.eps').*cosh(s.r').*sinh(s.r')).*(s.t').^2;
cp=1-sqrt(1./((n+1).^2-m.^2));

Ord = s.CO{s.k} - 1;
NM = s.M - Ord;

C = zeros(NM, 1);

for i = 1:NM
    C(i,:) = real(prod(cp(i:i+Ord,:),1));
end 
end                                         %End ncc function