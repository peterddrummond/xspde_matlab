function C = Delx(a,p)
% C = Delx(a,p) generates the x quadrature standard deviation for qsim
% follows conventions from the entanglement sction

C = 2 - p.phase;                                %ordering correction
C = C + (a(1:p.m,:) + a(p.m+1:2*p.m,:)).^2;       %x amplitude squared
C = sqrt(real(mean(C,2)));                       %Average over ensemble
end                                              %end x function