function d  =  xd1(a,varargin) 
%   a = XD1(a,dir,r) calculates a first spatial derivative of any field a
%   using second order (in dx) central finite difference methods.  
%   Output: first derivative  d.
%   Boundary values are specified in the field, except for periodic case
%   Derivatives at boundaries are for measurements, not for propagation
%   and are set equal to interior derivatives to give an approximate value
%   If given only two inputs, takes derivative in x direction
%   For three inputs, it takes transverse derivative in axis  'dir' = 2,3..
%   All boundary types and values can be set individually per field index.
%   They apply to any transverse dimension, and can change dynamically.
%   Boundary types are defined through r.boundaries{d}(i,b)
%   per dimension (d=1,2..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d is the transverse space dimension, ie, not counting time
%   (a) r.boundaries{d}(i,j)   = 0 gives the default, periodic 
%   (b) r.boundaries{d}(i,j) = -1 gives Neumann prescribed derivative
%   (c) r.boundaries{d}(i,j)= 1 gives Dirichlet, prescribed field
%   Boundary values are set on field i through boundvalue(i,..,j,..).
%   This array is dynamically returned by r.boundfun{dir}(a,r)
%   The index list of the boundary values matches the field index list.
%   Values are set sequentially by the dimension, and can be overwritten
%   First dimension of a is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License 
                                                            
if nargin < 2 || nargin >3
    error('Function xd1 takes two or three arguments, not %d', nargin)
end
sz = size(a);                            %%get original size to return
if nargin == 2
    r = varargin{1};                     %% r is the second argument
    dir = 2;                             %% x-derivative
else
    dir =varargin{1};                    %% direction is input
    r = varargin{2};                     %% r is the last argument
end
if dir < 2 || dir > r.dimension
    error('Derivative direction %d in xd1 is not supported',dir)
end
sz2 = prod(sz(2:end));
if r.nlattice ~= sz2
    error('xd1 expects lattice size of %d, not %d',r.nlattice,sz2)
end
shape(1) =  sz(1);
shape(2) =  prod(r.d.int(1:dir-1)); 
shape(3) =  r.d.int(dir); 
shape(4) =  prod(r.d.int(dir+1:end)); 
a = reshape(a,shape);                    %%Unflatten lattice
d = zeros(shape);
dx = r.dx(dir);
dx2 = 2*dx;
d(:,:,2:end-1,:) = (a(:,:,3:end,:)-a(:,:,1:end-2,:))/dx2;
shape(3) =  2;                           %%for boundary values
boundvalue=r.boundfun(a,dir,r);
boundvalue = reshape(boundvalue,shape);
for i=1:sz(1)
  if r.boundaries{dir}(i,1)  == 1        %%Dirichlet, prescribed field
    d(i,:,1,:) = d(i,:,2,:);
  elseif r.boundaries{dir}(i,1) == -1
        d(i,:,1,:) = boundvalue(i,:,1,:);%%Neuman, prescribed derivative
  else
    d(i,:,1,:) = (a(i,:,2,:)-a(i,:,end,:))/dx2;%%periodic case
  end
  if r.boundaries{dir}(i,2) == 1        %%Dirichlet, prescribed field
    d(i,:,end,:) = d(i,:,end-1,:);   
  elseif r.boundaries{dir}(i,2) == -1   %%Neuman, prescribed derivative
    d(i,:,end,:) = boundvalue(i,:,end,:);
  else
    d(i,:,end,:) = (a(i,:,1,:)-a(i,:,end-1,:))/dx2;%%periodic case
  end
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd1 function