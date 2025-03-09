function C = K1c(p)
% C = K1C(p); generates comparison click probabilities in one partition,
% assumed to be of length p.modes

n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
k = 1-sqrt(1./((n+1).^2-m.^2));
C = xchooseft1(p.modes, k);
C = reshape(C,[1,1,numel(C)]);

end 
