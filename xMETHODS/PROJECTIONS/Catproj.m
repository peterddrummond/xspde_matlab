function a  =  Catproj(d,a,n,r) 
%   a = Catproj(d,a,n,r) computes projections for a catenoid
%   The projection constraint function is:
%   f = (x_1)^2 + (x_2)^2 - (sinh(x_3))^2 - 1 = 0
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   xcatproj needs: 
%          r.iterations - defines how many normal iterations used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch n
case 0
  a =  2*[a(1,:);a(2,:);...
        -sinh(a(3,:)).*cosh(a(3,:))];           %%Returbs gradient vector
case 1
  gr =  2*[a(1,:);a(2,:);...
        -sinh(a(3,:)).*cosh(a(3,:))];           %%Calculate gradient vector
  a  = d-gr.*sum(d.*gr,1)./sum(gr.*gr,1);       %%Get projected derivative
case 2
  for i = 1:r.iterations
    gr =  2*[a(1,:);a(2,:);...
        -sinh(a(3,:)).*cosh(a(3,:))];           %%Calculate gradient vector  
    l = -constr(a,r)./sum(gr.*gr,1);            %%Get Lagrange multiplier
    a = a + l.*gr;                              %%Project the field again
  end
  case 3
  a =  constr(a,r);
otherwise
  error('Catproj expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,~)
%   a = constr(a,~) computes constraints for a catenoid
%   The projection constraint function is:
%   f = (x_1)^2 + (x_2)^2 - (sinh(x_3))^2 - 1 = 0
    f = a(1,:).^2+a(2,:).^2-sinh(a(3,:)).^2-1;  %%Calculate constraint
end      