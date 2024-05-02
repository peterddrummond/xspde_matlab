function d1  =  D1(a,varargin) 
%   d1 = D1(a[,d,j,c],p) calculates first spatial derivatives
%   using central finite differences for fields or averages  
%   Output: d1, the first derivative along dimension d. Returns:
%   (i)   two inputs:   all derivative indices in dimension 2
%   (ii)  three inputs: all derivative indices in dimension d
%   (iii) four inputs:  derivatives in dimension d and index list j
%   (iv) five inputs:   derivatives in dimension d, index list j, cell c.
%   If j is input as a vector, indices in the list are returned in order.
%   In this case, the output derivative may have a different size to a.
%   If p.indext = 0, there is no time index: used for derivative functions.
%   If p.indext = 1, there is a time index: used for functions of averages.
%   The field is assumed periodic except for each Neumann boundary
%   For Neumann boundaries, a cell index c is required unless c = 1
%   Requires p.boundval to be initialized before it is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: deriv, observe, output
%   Needs: p.dx, p.indext, p.boundaries, p.boundval
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sx = size(a);                                    % get original size
irange = 1:sx(1);                                % set range to all indices
p = varargin{end};                               % p is the last argument
d = varargin{1};                                 % dimension is input
c = 1;
switch nargin
  case 2
    d = 2;                                       % x-derivative
  case 3
  case 4
	irange = varargin{2};                        % component list is input
  case 5  
	irange = varargin{2};                        % component list is input
    c = varargin{3};                             % cell index is input
  otherwise
    error('Function D1 takes two to five arguments, not %d', nargin)
end
d = d + p.indext;                                % get dimension index
if p.indext == 0 && any(p.boundaries{c,d}(:) == -1)
  bval = p.boundval{c,d};
end
sx(1) = length(irange);
d1    = zeros(sx);                               % initialize output = zero
j     = 1;                                       % initial derivative index
for i = irange                                   % loop over field indices
  ind = p.ind;                                   % get index cell vector
  ind{1} = i;                                    % set first index cell
  ai  = a(ind{:});
  ind{1} = j;
  d1(ind{:}) = (circshift(ai,-1,d)-circshift(ai,1,d))/(2*p.dx(d));
  if p.indext == 0 
    if p.boundaries{c,d}(i,1) == -1              % Use Neumann at low end?
      ind{d} = 1;
      d1(ind{:}) = bval{i,1};
    end                                          % end if low end
    if p.boundaries{c,d}(i,2) == -1              % Use Neumann at high end?
      ind{d} = sx(d);
      d1(ind{:}) = bval{i,2};
    end                                          % end if high end
  end                                            % end if indext
  j = j + 1;                                     % increment deriv index
end                                              % end loop on field index
end                                              % end D1 function