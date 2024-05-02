function a  =  Quadprojc(d,a,n,r) 
%   a = Quadproj(d,a,n,r) computes projections for a complex quadratic.
%   The projection constraint function is, for hermitian qc:
%   f = \sum_i,j conj(a_i) * qc_{i,j} * (a_j) - 1 = 0
%   options available:
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   needs: r.qc as a constraint matrix, 
%          r.iterations defines how many normal iterations used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g  =  r.qc;                                   %%Quadratic constraint matrix
switch n
case 0
  a = 2.*g*a;                                 %%Calculate gradient vector 
case 1
  gr =  2.*g*a;                               %%Calculate gradient vector
  a  = d-gr.*sum(d.*conj(gr),1)./sum(gr.*conj(gr),1);%%Projected derivative
case 2
  for i = 1:r.iterations
      gr =  2.*g*a;                           %%Calculate gradient vector
      l = -constr(a,r)./sum(gr.*conj(gr),1);  %%Get Lagrange multiplier
      a = a+l.*gr;                            %%Project the field again
  end  
  case 3
  a =  constr(a,r);
otherwise
  error('expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,r)
    f = sum(conj(a).*(r.qc*a),1)-1;                 %%Calculate constraint
end           