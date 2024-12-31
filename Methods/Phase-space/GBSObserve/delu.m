function C = delu(a,p)
% C = delu(a,p) generates the u quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.m-1);
x = a(1:p.m,:) + a(p.m+1:2*p.m,:);
u = x(1,:) - sum(x(2:p.m,:),1)*fact;
C = sqrt(real(mean(u.^2,2)) -2*(p.phase-2));
end                                              %end du function