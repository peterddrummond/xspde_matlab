function C = Kmc(p)
% C = KMC(p); generates comparison O-th order click probability 

O = p.correl{p.noutput};                                    %list of moments
Ol = length(O);                                  %length of graph axis
n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
cp = 1-sqrt(1./((n+1).^2-m.^2));
C = zeros([1,1,Ol,1]);                                 %initialize correlation
for i = 1:Ol
    C(1,1,i,1) = prod(cp(1:O(i)));                   %Get O-th order click prob
end
end 
