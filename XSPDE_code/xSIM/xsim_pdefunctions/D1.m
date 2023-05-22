  function d1  =  D1(a,varargin) 
%   d1 = D1(a[,d,j],p) calculates first spatial derivatives
%   using central finite differences, for either fields or averages  
%   Output: d1, the first derivative along the axis dimension d.
%   (i)  two inputs: returns all component derivatives in x:  'd = 2'
%   (ii) three inputs: returns all derivatives in dimension 'd' = 2,3..
%   (ii) four inputs: returns a derivative just for components in list 'j'.
%   If j is a vector, indices in the list are returned sequentially.
%   In this case, the output derivative may have a different size to a.
%   Boundary values are used as specified in the parameter struucture p
%   All boundary types and values can be set individually per component.
%   They apply to any transverse dimension, and can change dynamically.
%%%%%%%%%%%%%%%%%%%%
%   BOUNDARY TYPES %
%%%%%%%%%%%%%%%%%%%%
%       These are defined by the input parameter p.boundaries{d}(j,b)
%   per dimension (d = 2,..), index (j = 1,2..) and boundary b=(1,2).
%   Here b = 1 is the lower boundary index, b = 2 is the upper index. 
%   As elsewhere d is the dimension index of the specified boundary type. 
%   No boundary is set for d = 1, as D1 only sets spatial boundaries. 
%   (a) p.boundaries{d}(j,b)  = 0  is the default, periodic 
%   (b) p.boundaries{d}(j,b)  = -1 is Neumann/Robin, prescribed derivative
%   (c) p.boundaries{d}(j,b)  = 1  is Dirichlet, prescribed field
%%%%%%%%%%%%%%%%%%%%%
%   BOUNDARY VALUES %
%%%%%%%%%%%%%%%%%%%%%
%   These are set dynamically on field j using an array bval(j,..,b,..).
%   This array is returned by the user function p.boundfun(a,d,p)
%   The indices of the returned boundary values matches the field indices,
%   except that there only two index values are needed in the direction d.
%   Values are set sequentially by the dimension index d = 2,3,..
%   If p.indext = 0, there is no time index and boundary values are used.
%   If p.indext = 1 for observed averages, boundary values are set to zero.
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License 
                                                            
irange = 0;
switch nargin
    case 2
        dir = 2;                                 %% x-derivative
    	p = varargin{1};                         %% p is the last argument
    case 3
    	dir = varargin{1};                       %% direction is input
    	p = varargin{2};                         %% p is the last argument
    case 4
        dir = varargin{1};                       %% direction is input
	    irange = varargin{2};                    %% component list is input
        p = varargin{3};                         %% p is the last argument
    otherwise
        error('Function D1 takes two to four arguments, not %d', nargin)
end
sz = size(a);                                    %% original size to return
in = p.indext;                                   %% get time-index switch
s(1) =  sz(1);
s(2) =  prod(sz(2:in+dir-1)); 
s(3) =  sz(in+dir); 
s(4) =  prod(sz(in+dir+1:end));
if p.indext == 0                                 %%propagating, use b.c.
    bval = p.boundfun(a,dir,p);
    bval = reshape(bval,[s(1:2),2,s(4)]);
else                                             %%post-processing, no b.c.
    bval = zeros([s(1:2),2,s(4)]);
end
a = reshape(a,s);                                %%reshape lattice
if isequal(irange,0)                             %%no input component list
    irange = 1:sz(1);                            %%set range to components
end
sz(1) = length(irange);
s(1)  = length(irange);
d1 = zeros(s);
dx = p.dx(dir);
dx2 = 2*dx;
j1 = 1;
for j=irange
  d1(j1,:,2:end-1,:) = (a(j,:,3:end,:)-a(j,:,1:end-2,:))/dx2;
  if p.boundaries{dir}(j,1)  == 1                %%Dirichlet, field value
    d1(j1,:,1,:) = (a(j,:,2,:) - bval(j,:,1,:))/dx2;
  elseif p.boundaries{dir}(j,1) == -1
    d1(j1,:,1,:) = bval(j,:,1,:);                %%Neuman/Robin, derivative
  else
    d1(j1,:,1,:) = (a(j,:,2,:)-a(j,:,end,:))/dx2;%%periodic case
  end
  if p.boundaries{dir}(j,2) == 1                 %%Dirichlet, field value
    d1(j1,:,end,:) = (bval(j,:,end,:)-a(j,:,end-1,:))/dx2;   
  elseif p.boundaries{dir}(j,2) == -1            %%Neuman/Robin, derivative
    d1(j1,:,end,:) = bval(j,:,end,:);
  else
    d1(j1,:,end,:) = (a(j,:,1,:)-a(j,:,end-1,:))/dx2; %%periodic case
  end
  j1=j1+1;
end
d1 = reshape(d1,sz);                       %%reshape to matrix
end                                      %%end xd1 function