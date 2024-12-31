function C = nc(p)
% C = NC(p); generates comparison numbers of  clicks per channel
% assumes output unchanged from input, apart from transmission factor t
% eg, identity matrix or uniform thermal+unitary

C = (sinh(p.sqz).^2+p.alpha.^2).*p.tr.^2;
C = reshape(C,[1,1,numel(C)]);
end                                         %End nc function