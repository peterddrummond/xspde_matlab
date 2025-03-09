function C = K1(a,p)
% C = K1(a,p); generates click probabilities in a single partition

if p.phase ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end
np2 = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);         %normal ordered numbers
cp  = 1-exp(-np2);                          %click projector
C = xchooseft1(p.modes,cp);
sz = size(C);
C = reshape(C, [sz(2), p.ensembles(1)]);
end                                         %End cprob1 function