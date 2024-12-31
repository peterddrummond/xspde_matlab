function C = delv(a,p)
% C = delv(a,p) generates the v quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.m-1);
y = (a(1:p.m,:) - a(p.m+1:2*p.m,:))/1i;
v = y(1,:) + sum(y(2:p.m,:),1)*fact;
C = sqrt(real(mean(v.^2,2)) -2*(p.phase-2));
end                                              %end du function