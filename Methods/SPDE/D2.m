function d2  =  D2(a,varargin) 
%   d2 = D2(a[,d,c,j],p) calculates second spatial derivatives of input a
%   using central finite differences for either fields or averages  
%   Output: d2, the second derivative along dimension d. Returns:
%   (i)   two inputs:   all derivatives in dimension 2, cell 1
%   (ii)  three inputs: all derivatives in dimension d, cell 1
%   (iii) four inputs:  all derivatives in dimension d, cell c
%   (iv)  five inputs:  derivatives in dimension d, cell c, indices j.
%   If j is input as a vector, indices in the list are returned in order.
%   In this case, the output derivative may have a different size to a.
%   If p.indext = 0, there is no time index: used for derivative functions.
%   If p.indext = 1 or c < 1, no boundary values used; D2 is continuous
%   For specified boundaries, c is needed unless c = 1 or indext = 1
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
    c = varargin{2};                             % cell index is input
  case 5  
	irange = varargin{3};                        % component list is input
    c = varargin{2};                             % cell index is input
  otherwise
    error('Function D2 takes two to five arguments, not %d', nargin)
end
dx = p.dx(d);                                    % space step  
d = d + p.indext;                                % get dimension index
if p.indext == 0 && c > 0 
  if any(p.boundaries{c,d}(:) == -1)
    bval = p.boundval{c,d};
  end
else
  c = 0;
end
sx(1) = length(irange);
d2    = zeros(sx);                               % initialize output = zero
j     = 1;                                       % initial derivative index
for i = irange                                   % loop over field indices
  in = p.ind; out = in;                          % initialize index cells
  in{1} = i;                                     % set first index to i
  ai  = a(in{:});                                % get input slice at i
  out{1} = j;                                    % set output slice index 
  d2(out{:}) = (circshift(ai,-1,d)+circshift(ai,1,d)-2*ai)/dx^2;
  if c < 1 || p.boundaries{c,d}(i,1) ~= 0        % if not periodic   
    in{d}  = 2; out{d} = 1; in{1} = 1;           % set low end indices
    if c > 0 && p.boundaries{c,d}(i,1) == -1     % Neumann low end
      a2 = ai(in{:});                            % get a2 with in{d} = 2
      in{d} = 1;                                 % reset in{d} = 1
      d2(out{:}) = 2*((a2-ai(in{:}))/dx-bval{i,1})/dx; %get boundary term 1
    else                                         % not Neumann or c = 0
      out2 = out; out2{d} = 2;                   % set index d to 2 for d2
      d2(out{:}) = d2(out2{:});                  % make a continuous D2
    end                                          % end if not Neuman
    in{d}  = sx(d)-1; out{d} = sx(d);            % high end indices
    if c > 0 && p.boundaries{c,d}(i,2) == -1     % Neumann high end?
      a2 = ai(in{:});                            % get a2; in{d} = sx(d)-1
      in{d} = sx(d);                             % reset in{d} = sx(d)
      d2(out{:}) = 2*(bval{i,2}+(a2-ai(in{:}))/dx)/dx; %get boundary term 2
    else                                         % not Neumann or c = 0
      out2 = out; out2{d} = sx(d)-1;             % set index d to 2 for d2
      d2(out{:}) = d2(out2{:});                  % make a continuous D2
    end                                          % end Neumann high end?
  end                                            % end if not periodic
  j = j + 1;                                     % increment deriv index
end                                              % end loop on field index
end                                              % end D2 function