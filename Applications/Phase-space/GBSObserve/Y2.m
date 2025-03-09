function C = Y2(a,p)
% C = Y2(a,p) generates the y quadrature second moment for qsim

C = 2 - p.phase;                                 %ordering correction
if p.phase == 1 && size(a,1) == 2*p.modes        % if positive P
  C = C - (a(1:p.modes,:)-a(p.modes+1:2*p.modes,:)).^2; % +P case
else                                             % else not +P
  C = C - (a-conj(a)).^2;                        % add ordering correction
end                                              % end if normal-ordered 
end                                              %end x function