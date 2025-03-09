function a  =  Polproj(d,ac,n,varargin) 
%   a = Polproj(d,ac,n,(c,)p) computes projections for a polynomial
%   The projection constraint function is:
%   f = \sum_j vc_j * (a_j)^orderpol - 1 = 0
%   options available:
%   n = 0 - returns a tangential vector at location a
%   n = 1 - returns a tangential projection of d at location a
%   n = 2 - returns a normal projection of a to the constraint manifold
%   n = 3 - returns the constraint function
%   Here c is optional: if present, ac is a cell array and c the cell index
%   Needs: 
%          p.orderpol defines polynomial order
%          p.vc defines polynomial coefficient vector
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
g  =  p.vc{c};                                   % coefficient column vector
o  =  p.orderpol{c};
switch n
case 0
  if isempty(g)
    a = 0;
  else
    a = o*g.*a.^(o-1);                             % Calculate gradient vector 
  end
case 1
  if isempty(g)
    a = d;
  else
    gr = o*g.*a.^(o-1);                            % Calculate gradient vector 
    a  = d - gr.*sum(d.*gr,1)./sum(gr.*gr,1);      % Get projected derivative
  end
case 2
  if ~isempty(g)
    for i = 1:p.iterations
      gr =  o*g.*a.^(o-1);                       % Calculate gradient vector 
      l = -constr(a,c,p)./sum(gr.*gr,1);         % Get Lagrange multiplier
      a = a+l.*gr;                               % Project the field again
    end
  end  
case 3
  a =  constr(a,c,p);
otherwise
  error('expects n = 0,1,2 or 3, not %d',n);
end
end

function f  =  constr(a,c,p)
if isempty(p.vc{c})
    f = 0;
else
    f = sum((p.vc{c}).*a.^p.orderpol{c},1)-1;    % Calculate constraint
end
end           