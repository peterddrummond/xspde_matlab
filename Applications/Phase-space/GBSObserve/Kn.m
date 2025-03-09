function C = Kn(a,p)
% C = KN(a,p); generates click probabilities for an n-fold partition
% The partition is defined by the quantity p.part{p.nobserve}
% This is EITHER  a cell array of index vectors
% OR, a vector of partition lengths

if p.phase ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end
n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);    %normal ordered numbers
k  = 1-exp(-n);                                  %click projector
if iscell(p.part{p.nobserve})                    %check if cell partition 
    C = xchooseftnc(p.part{p.nobserve},k);       %get click probabilities
else                                             %not cell partition 
    C = xchooseftn(p.part{p.nobserve},k);        %get click probabilities
end                                              %check if cell partition 
end                                              %end kn