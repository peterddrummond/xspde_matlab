function a  =  Quadproj(d,ac,n,varargin) 
%   a = Quadproj(d,ac,n,(c,)p) computes projections for a quadratic.
%   The projection constraint function is:
%   f = \sum_i,j a_i * qc_{i,j} * (a_j) - 1 = 0
%   options available, if p.qc{c} is not empty:
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   Here c is optional: if present, ac is a cell array and c the cell index
%   needs: p.qc as a cell array of constraint matrices, 
%          p.iterations defines how many normal iterations used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 4 
    a = ac;
    c = 1;
    p = varargin{1};
else
    c = varargin{1};
    p = varargin{2};
    a = ac{c};
end
g  =  p.qc{c};                                %%Quadratic constraint matrix
switch n
case 0
  if isempty(g)
    a = 0;
  else
    a = 2.*g*a;                                 %%Calculate gradient vector
  end 
case 1
  if isempty(g)
    a = d;
  else
    gr =  2.*g*a;                               %%Calculate gradient vector
    a  = d-gr.*sum(d.*gr,1)./sum(gr.*gr,1);     %%Get projected derivative
  end
case 2
  if ~isempty(g)
    for i = 1:p.iterations
      gr =  2.*g*a;                           %%Calculate gradient vector
      l = -constr(a,c,p)./sum(gr.*gr,1);        %%Get Lagrange multiplier
      a = a+l.*gr;                            %%Project the field again
    end
  end  
case 3
  a =  constr(a,c,p);
otherwise
  error('expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,c,p)
  if isempty(p.qc{c})
    f = 0;
  else
    f = sum(a.*(p.qc{c}*a),1)-1;                 %%Calculate constraint
  end
end           