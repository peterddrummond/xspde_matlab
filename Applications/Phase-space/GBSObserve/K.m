function C = K(a,p)
% C = K(a,p); generates mean numbers of clicks per channel for qsim

if p.phase ~= 1
    error('xSPDE supports click probabilities only for +P phasespace = 1');
end
n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);        %normal ordered numbers
C  = 1-exp(-n);                                  %click projector
end                                              %End nc function