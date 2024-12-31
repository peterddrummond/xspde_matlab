function C = n1thc(p)
% C = N1THC(p); generates comparison count probabilities in one partition
% of length p.modes, with thermal inputs,  p.thermal=1, and unitary transmission
if p.thermal ~= 1
    error('n1c uses thermal noise inputs');
end
p1  = 1/(1+mean((sinh(p.sqz).*p.tr).^2,2));
C = p.modes*log(p1)+zeros(p.max+1,1);
lp2 = log(1-p1);
for m = 1:p.max
    C(m+1,1) = C(m,1)+log(1+(p.M-1)/m)+lp2;
end
C = exp(C);
