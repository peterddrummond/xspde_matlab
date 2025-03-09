  function d1  =  D1(a,varargin) 
%   d1 = D1(a[,d,c,j],p) calculates first spatial derivatives
%   using central finite differences, for fields or averages  
%   Output: d1, the first derivative along the axis dimension d. Returns if
%   (i)  two inputs:  all derivatives in x:  'd = 2', assumes cell(1)
%   (ii) three inputs: all derivatives in dimension 'd', assumes cell(1) 
%   (iii) four inputs: derivatives of all components of cell(c).
%   (iv) five inputs: derivatives of components in list 'j' of cell (c).
%   If j is a vector, indices in the list are returned sequentially.
%   In this case, the output derivative may have a different size to a.
%   Boundary values are used as specified in the parameter structure p
%   All boundary types and values can be set individually per component.
%   They apply to any transverse dimension, and can change dynamically.
%   The cell index is used to define the cell boundary value that is set

%%%%%%%%%%%%%%%%%%%%
%   BOUNDARY TYPES %
%%%%%%%%%%%%%%%%%%%%
%
%       These are defined by the input parameter p.boundaries{c,d}(j,b)
%   per cell (c), dimension (d = 2,.), index (j = 1,..), boundary b=(1,2).
%   Here b = 1 is the lower boundary index, b = 2 is the upper index. 
%   As elsewhere d is the dimension index of the specified boundary type. 
%   No boundary is set for d = 1, as D1 only sets spatial boundaries. 
%   (a) p.boundaries{c,d}(j,b) = 0  is the default, periodic 
%   (b) p.boundaries{c,d}(j,b) = -1 is Neumann/Robin, prescribed derivative
%   (c) p.boundaries{c,d}(j,b) = 1  is Dirichlet, prescribed field

%%%%%%%%%%%%%%%%%%%%%
%   BOUNDARY VALUES %
%%%%%%%%%%%%%%%%%%%%%
%
%   These are set using an array  p.boundval{c,d}{j,b}.
%   The indices of the returned boundary values matches the field indices,
%   except that there only one index value is needed in the direction d.
%   Values are set sequentially by the dimension index d = 2,3,..
%   If p.indext = 0, there is no time index and boundary values are used.
%   If p.indext = 1 for observed averages, boundary values are set to zero.
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License 
                                                            
irange = 0;
c=1;
switch nargin
    case 2
        d = 2;                                   %% x-derivative
    	p = varargin{1};                         %% p is the last argument
    case 3
    	d = varargin{1};                         %% direction is input
    	p = varargin{2};                         %% p is the last argument
    case 4
        d = varargin{1};                         %% direction is input
	    c = varargin{2};                         %% cell is input
        p = varargin{3};                         %% p is the last argument
    case 5
        d = varargin{1};                         %% direction is input
        c = varargin{2};                         %% cell is input
	    irange = varargin{3};                    %% component list is input
        p = varargin{4};                         %% p is the last argument
    otherwise
        error('Function D1 takes two to five arguments, not %d', nargin)
end
sz = size(a);                                    %% original size to return
in = p.indext;                                   %% get time-index switch
s(1) =  sz(1);
s(2) =  prod(sz(2:in+d-1)); 
s(3) =  sz(in+d); 
s(4) =  prod(sz(in+d+1:end));
a = reshape(a,s);                                %%reshape lattice
if isequal(irange,0)                             %%no input component list
    irange = 1:sz(1);                            %%set range to components
end
sz(1) = length(irange);
sb = s;
sb(3) = 1;
sb(1) = 1;
s(1)  = sz(1);
if p.indext == 0 && c > 0                        %%propagating, use b.c.
    if any(p.boundaries{c,d}(:) == -1)
      bval = p.boundval{c,d};
    end
else                                             %%post-processing, no b.c.
    c = 0;
end
d1 = zeros(s);
dx = p.dx(d);
dx2 = 2*dx;
j1 = 1;
for j=irange
  aj = a(j,:,:,:);
  d1(j1,:,:,:) = (circshift(aj,-1,3)-circshift(aj,1,3))/dx2; % periodic
  if c == 0                                      %%No boundaries specified
    d1(j1,:,1,:) = d1(j1,:,2,:);
    d1(j1,:,end,:) = d1(j1,:,end-1,:);
  else
    if p.boundaries{c,d}(j,1)  == - 1
      d1(j1,:,1,:) = reshape(bval{j,1},sb);      %%Neuman/Robin, derivative
    end
    if p.boundaries{c,d}(j,2)  == - 1         
      d1(j1,:,end,:) = reshape(bval{j,2},sb);    %%Neuman/Robin, derivative
    end
  end
  j1=j1+1;
end
d1 = reshape(d1,sz);                             %%reshape to input size
end                                              %%end D1 function