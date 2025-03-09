function C = Dely(a,p)
% C = Dely(a,p) generates the y quadrature standard deviation

C = 2 - p.phase;                                %ordering correction
C = C - (a(1:p.m,:) - a(p.m+1:2*p.m,:)).^2;      %Y amplitude squared
C = sqrt(real(mean(C,2)));                       %Average over ensemble
end                                              %end x function