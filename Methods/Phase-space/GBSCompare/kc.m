function C = kc(p)
% C = KC(p); generates comparison numbers of clicks per channel
% assumes output unchanged from input. 
% eg, identity matrix or uniform thermal+unitary
n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
C = 1-sqrt(1./((n+1).^2-m.^2));
C = reshape(C,[1,1,numel(C)]);
end                                         %End ncc function