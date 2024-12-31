function C = delu2(a,p)
% C = delu(a,p) generates the u quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.modes-1);
x = a(1:p.modes,:) + a(p.modes+1:2*p.modes,:);
u = x(1,:) - sum(x(2:p.modes,:),1)*fact;
C = real(u.^2) -2*(p.phase-2);
end                                              %end du function