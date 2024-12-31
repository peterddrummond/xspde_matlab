function C = deltot(a,p)
% C = delv(a,p) generates the v quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.m-1);
p = (a(1:p.m,:) - a(p.m+1:2*p.m,:))/1i;
v = p(1,:) + sum(p(2:s.m,:),1)*fact;
C = sqrt(real(v.^2) -2*(p.phase-2));
end                                              %end du function