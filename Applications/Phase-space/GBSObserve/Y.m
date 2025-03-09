function C = Y(a,p)
% C = Y(a,p) generates the y quadrature 

C = -i*(a(1:p.modes,:)+ a(p.modes+1:2*p.modes,:));    %y amplitude
end                                                   %end function