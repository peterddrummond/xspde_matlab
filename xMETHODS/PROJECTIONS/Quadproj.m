function a  =  Quadproj(d,ac,n,varargin) 
%   a = Quadproj(d,ac,n,(c,)p) computes projections for a quadratic.
%   The projection constraint function is:
%   f = \sum_i,j a_i * qc_{i,j} * (a_j) - 1 = 0
%   options available:
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   Here c is optional: if present, ac is a cell array and c the cell index
%   needs: p.qc as a constraint matrix, 
%          p.iterations defines how many normal iterations used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 4 
    a = ac;
    p = varargin{1};
else
    c = varargin{1};
    p = varargin{2};
    a = ac{c};
end
g  =  p.qc;                                   %%Quadratic constraint matrix
switch n
case 0
  a = 2.*g*a;                                 %%Calculate gradient vector 
case 1
  gr =  2.*g*a;                               %%Calculate gradient vector
  a  = d-gr.*sum(d.*gr,1)./sum(gr.*gr,1);     %%Get projected derivative
case 2
  for i = 1:p.iterations
      gr =  2.*g*a;                           %%Calculate gradient vector
      l = -constr(a,p)./sum(gr.*gr,1);        %%Get Lagrange multiplier
      a = a+l.*gr;                            %%Project the field again
  end  
  case 3
  a =  constr(a,p);
otherwise
  error('expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,p)
    f = sum(a.*(p.qc*a),1)-1;                 %%Calculate constraint
end           