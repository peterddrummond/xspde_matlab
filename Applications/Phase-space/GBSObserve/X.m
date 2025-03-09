function C = X(a,p)
% C = X(a,p) generates the x quadrature 

C = (a(1:p.modes,:)+ a(p.modes+1:2*p.modes,:));           %y amplitude
end                                              %end x function