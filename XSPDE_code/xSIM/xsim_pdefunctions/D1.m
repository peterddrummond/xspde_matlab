  function d  =  D1(a,varargin) 
%   a = D1(a[,dir, component],r) calculates first spatial derivatives
%   using central finite differences.  
%   Output: d, the first derivative of along the axis dimension dir.
%   Boundary values are specified in the field, except in the periodic case
%   If given only two inputs, takes derivatives in x direction:  'dir=2'
%   For three inputs, it takes transverse derivatives in axis  'dir' = 2,3..
%   For four inputs, it returns a single derivative for 'i=1,2,3..'.
%   All boundary types and values can be set individually per component.
%   They apply to any transverse dimension, and can change dynamically.
%   Boundary types are defined through r.boundaries{d}(i,b)
%   per dimension (d = 2,..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d is the dimension index of the specified boundary type. No
%   boundary is set for d = 1, as xd1 only sets spatial boundaries. 
%   (a) r.boundaries{d}(i,j)  = 0 gives the default, periodic 
%   (b) r.boundaries{d}(i,j)  = -1 gives Neumann, prescribed derivative
%   (c) r.boundaries{d}(i,j)  = 1 gives Dirichlet, prescribed field
%   Boundary values are set on field i through boundvalue(i,..,j,..).
%   This array is dynamically returned by r.boundfun{dir}(a,r)
%   The index list of the boundary values matches the field index list.
%   Values are set sequentially by the dimension, and can be overwritten
%   First dimension of a is the field index, last dimension is the ensemble
%   If r.indext = 0, there is no time index, otherwise there is one.
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License 
                                                            
comp = 0;
switch nargin
    case 2
    	r = varargin{1};                     %% r is the second argument
    	dir = 2;                             %% x-derivative
    case 3
    	dir =varargin{1};                    %% direction is input
    	r = varargin{2};                     %% r is the last argument
    case 4
        dir =varargin{1};                    %% direction is input
	    comp = varargin{2};                  %% component is input
        r = varargin{3};                     %% r is the last argument
    otherwise
        error('Function D1 takes two to four arguments, not %d', nargin)
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
a = reshape(a,s);                             %%flatten lattice
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
d = zeros(s);
dx = r.dx(dir);
dx2 = 2*dx;
i1 = 1;
for i=irange
  d(i1,:,2:end-1,:) = (a(i,:,3:end,:)-a(i,:,1:end-2,:))/dx2;
  if r.boundaries{dir}(i,1)  == 1        %%Dirichlet, prescribed field
    d(i1,:,1,:) = (a(:,:,2,:) - bval(i,:,1,:))/dx2;
  elseif r.boundaries{dir}(i,1) == -1
    d(i1,:,1,:) = bval(i,:,1,:);    %%Neuman, prescribed derivative
  else
    d(i1,:,1,:) = (a(i,:,2,:)-a(i,:,end,:))/dx2;%%periodic case
  end
  if r.boundaries{dir}(i,2) == 1        %%Dirichlet, prescribed field
    d(i1,:,end,:) = (bval(i,:,end,:)-a(i,:,end-1,:))/dx2;   
  elseif r.boundaries{dir}(i,2) == -1   %%Neuman, prescribed derivative
    d(i1,:,end,:) = bval(i,:,end,:);
  else
    d(i1,:,end,:) = (a(i,:,1,:)-a(i,:,end-1,:))/dx2; %%periodic case
  end
  i1=i1+1;
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd1 function