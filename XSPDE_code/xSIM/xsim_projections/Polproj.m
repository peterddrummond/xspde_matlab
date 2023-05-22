function a  =  Polproj(d,a,n,r) 
%   a = Polproj(d,a,n,r) computes projections for a polynomial
%   The projection constraint function is:
%   f = \sum_j vc_j * (a_j)^orderpol - 1 = 0
%   options available:
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   needs: 
%          r.orderpol defines polynomial order
%          r.vc defines polynomial coefficient vector
%          r.iterproj defines how many normal iterations used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g  =  r.vc;                                   %%coefficient column vector
switch n
case 0
  a = r.orderpol*g.*a.^(r.orderpol-1);        %%Calculate gradient vector 
case 1
  gr = r.orderpol*g.*a.^(r.orderpol-1);       %%Calculate gradient vector 
  a  = d - gr.*sum(d.*gr,1)./sum(gr.*gr,1);   %%Get projected derivative
case 2
  for i = 1:r.iterproj
      gr =  r.orderpol*g.*a.^(r.orderpol-1);  %%Calculate gradient vector 
      l = -constr(a,r)./sum(gr.*gr,1);        %%Get Lagrange multiplier
      a = a+l.*gr;                            %%Project the field again
  end  
  case 3
  a =  constr(a,r);
otherwise
  error('expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,r)
    f = sum((r.vc).*a.^r.orderpol,1)-1;      %%Calculate constraint
end           