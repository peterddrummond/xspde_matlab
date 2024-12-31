function C = knc(p)
% C = KNC(p); generates comparison click probabilities 
% for an n-fold partition. 
n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
k = 1-sqrt(1./((n+1).^2-m.^2));
if iscell(p.part{p.noutput})                     %check if cell partition 
    C = xchooseftnc(p.part{p.noutput},k);              %get click probabilities
else                                             %not cell partition 
    C = xchooseftn(p.part{p.noutput},k);        %get click probabilities
end                                              %check if cell partition
C = reshape(C,[1,1,size(C)]);
end 
