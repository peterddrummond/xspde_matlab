function d2  =  D2(a,varargin) 
%   d2 = D2(a[,d,j],p) calculates second spatial derivatives of input a
%   using central finite differences for either fields or averages  
%   Output: d2, the second derivative along dimension d. Returns:
%   (i)   two inputs:   all derivatives in x dimension, (d = 2)
%   (ii)  three inputs: all derivatives in dimension d
%   (iii) four inputs:  derivatives in dimension d and index list j
%   (iv)  five inputs:  derivatives in dimension d, index list j, cell c.
%   If j is input as a vector, indices in the list are returned in order.
%   In this case, the output derivative may have a different size to a.
%   If p.indext = 0, there is no time index: used for derivative functions.
%   If p.indext = 1, there IS a time index: used for functions of averages.
%   The field is assumed periodic except for each Neumann boundary
%   For Neumann boundaries, c is needed unless c = 1 or indext = 1
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
c = 1;                                           % default cell
switch nargin
  case 2
    d = 2;                                       % x-derivative is default
  case 3                                         % no changes needed
  case 4
	irange = varargin{2};                        % component list is input
  case 5  
	irange = varargin{2};                        % component list is input
    c = varargin{3};                             % cell index is input
  otherwise
    error('Function D2 takes two to five arguments, not %d', nargin)
end
dx = p.dx(d);                                    % space step  
d = d + p.indext;                                % get dimension index
if p.indext == 0 && any(p.boundaries{c,d}(:) == -1)
  bval = p.boundval{c,d};
end
sx(1) = length(irange);
d2    = zeros(sx);                               % initialize output = zero
j     = 1;                                       % initial derivative index
for i = irange                                   % loop over field indices
  in = p.ind; inj = in;
  in{1} = i;                                % set first index to i
  ai  = a(in{:});                                % choose a slice at i
  in{1} = 1;                                     % reset first index of ai
  inj{1} = j;                               % index for a slice at j
  d2(inj{:}) = (circshift(ai,-1,d)+circshift(ai,1,d)-2*ai)/dx^2;
  if p.indext == 0 
    if p.boundaries{c,d}(i,1) == -1              % Use Neumann at low end?
      in{d}  = 2;                                % set index d to 2 for ai
      inj{d} = 1;                                % set index d to 1 for d2
      a2 = ai(in{:});                            % get ai with in{d} = 2
      in{d} = 1;                                 % set in{d} = 1
      d2(inj{:}) = 2*((a2-ai(in{:}))/dx-bval{i,1})/dx; %get boundary term 1
    end                                          % end if low end
    if p.boundaries{c,d}(i,2) == -1              % Use Neumann at high end?
      in{d}  = sx(d);inj{d} = sx(d);
      a2 = ai(in{:});
      in{d} = sx(d)-1;
      d2(inj{:}) = 2*(bval{i,2}-(a2-ai(in{:}))/dx)/dx;  %get boundary term 2
    end                                          % end if high end
  end                                            % end if indext
  j = j + 1;                                     % increment deriv index
end                                              % end loop on field index
end                                              % end D2 function