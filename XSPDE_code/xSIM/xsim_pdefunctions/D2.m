function d2  =  D2(a,varargin) 
%   d2 = D2(a[,d,j],p) calculates first spatial derivatives
%   using central finite differences, for either fields or averages  
%   Output: d2, the second derivative along the axis dimension d.
%   (i)  two inputs: returns all derivatives in x:  'd = 2'
%   (ii) three inputs: returns all derivatives in dimension 'd' = 2,3..
%   (ii) four inputs: returns a derivative for all components in list 'j'.
%   Note that if j is a vector, all indices in the list are returned.
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

comp = 0;
switch nargin
    case 2
    	r = varargin{1};                 %% r is the second argument
    	dir = 2;                         %% x-derivative
    case 3
    	dir =varargin{1};                %% direction is input
    	r = varargin{2};                 %% r is the last argument
    case 4
        dir =varargin{1};                %% direction is input
	    comp = varargin{2};              %% component is input
        r = varargin{3};                 %% r is the last argument
    otherwise
        error('Function D2 takes two to four arguments, not %d', nargin)
end
sz = size(a);                            %%get original size to return
in = r.indext;                           %% get time-index switch
s(1) =  sz(1);
s(2) =  prod(sz(2:in+dir-1)); 
s(3) =  sz(in+dir); 
s(4) =  prod(sz(in+dir+1:end));
if r.indext == 0 
    bval = r.boundfun(a,dir,r);
    bval = reshape(bval,[s(1:2),2,s(4)]);
else
    bval = zeros([s(1:2),2,s(4)]);
end
a = reshape(a,s);                        %%flatten lattice
irange = 1:sz(1);
if comp > 0 
    bval = bval(comp,:,:,:);
    irange = comp;
    if sz(1) > 1 
       a = a(comp,:,:,:);
       sz(1) = 1;
       s(1) = 1;
    end
end
d2 = zeros(s);                           %% Initialize derivative to zero
dx = r.dx(dir);                          %% Set the space step value
dx2 = dx^2;                              %% Set the space step squared
e = s(3);                                %% Set the extent
i1 = 1;                                  %% Initialize derivative index
for i=irange                             %% Loop over field indices
    d2(i1,:,2:e-1,:) = (a(i,:,3:e,:)+a(i,:,1:e-2,:)-2.*a(i,:,2:e-1,:))/dx2;
  if r.boundaries{dir}(i,1) == 1         %% Use Dirichlet at low end?
    d2(i1,:,1,:) = (a(i,:,2,:)-bval(i,:,1,:))/dx2;
  elseif r.boundaries{dir}(i,1) == -1    %% Use Neumann at low end?
    d2(i1,:,1,:) = 2*((a(i,:,2,:)-a(i,:,1,:))/dx-bval(i,:,1,:))/dx;
  else                                   %% Derivative is periodic
    d2(i1,:,1,:) = (a(i,:,2,:)-2*a(i,:,1,:)+a(i,:,e,:))/dx2;
  end                                    %% End lower boundary case
  if r.boundaries{dir}(i,2)  == 1        %% Use Dirichlet at high end?
    d2(i1,:,e,:) = (a(i,:,e-1,:)-bval(i,:,2,:))/dx2;
  elseif  r.boundaries{dir}(i,2)  == -1  %% Use Neumann at high end?
    d2(i1,:,e,:) = 2*(bval(i,:,2,:)-(a(i,:,e,:)-a(i,:,e-1,:))/dx)/dx;
  else                                   %% Derivative is periodic
    d2(i1,:,e,:) = (a(i,:,1,:)-2*a(i,:,e,:)+a(i,:,e-1,:))/dx2;
  end                                    %% End upper boundary case
i1 = i1+1;                               %% Increment derivative index
end                                      %% End loop over field indices
d2 = reshape(d2,sz);                     %% reshape to output matrix size
end                                      %% end D2 function