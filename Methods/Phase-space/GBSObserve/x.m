function C = x(a,p)
% C = X(a,p) generates the x quadrature for qsim

C = (a(1:p.modes,:)+ a(p.modes+1:2*p.modes,:));           %y amplitude
end                                              %end x function