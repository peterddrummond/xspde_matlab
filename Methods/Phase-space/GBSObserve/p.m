function C = p(a,p)
% C = P(a,p) generates the p quadrature for qsim

C = -i*(a(1:p.modes,:)+ a(p.modes+1:2*p.modes,:));           %y amplitude
end                                              %end x function