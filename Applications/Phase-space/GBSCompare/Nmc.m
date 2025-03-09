function C = Nmc(p)
% C = NMC(p); generates comparison O-th order moments of photon number

O = p.correl{p.noutput};                                    %list of moments
Ol = length(O);
n = (sinh(p.sqz').*p.tr').^2;
for i = 1:Ol
    C(1,1,i,1) = prod(n(1:O(i)));                    %Get O-th order moment
end
end 