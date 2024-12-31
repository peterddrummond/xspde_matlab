function C = delv2(a,p)
% C = delv2(a,p) generates the v quadrature variance  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.modes-1);
y = (a(1:p.modes,:) - a(p.modes+1:2*p.modes,:))/1i;
v = y(1,:) + sum(y(2:p.modes,:),1)*fact;
C = real(v.^2) -2*(p.phase-2);
end                                              %end du function