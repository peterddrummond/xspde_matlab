function d2  =  D2(a,varargin) 
%   d2 = D2(o,dir,r) calculates a second spatial derivative of any field a  
%   using second order (in dx) central finite difference methods.  
%   Output: second derivative  d2.
%   If two inputs, takes derivative in x direction
%   If three inputs, takes transverse derivative in axis  'dir' = 2,3..
%   Boundary conditions are in the field, except for periodic case
%   All boundary types and values can be set individually per field index.
%   They apply to any transverse dimension, and can change dynamically.
%   Boundary types are defined through r.boundaries{d}(i,b)
%   per dimension (d = 2,..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d is the dimension index of the specified boundary type. No
%   boundary can be set for d=1, as xd2 only sets spatial boundaries. 
%   (a) r.boundaries{d}(i,j)  = 0 gives the default, periodic 
%   (b) r.boundaries{d}(i,j)  = -1 gives Neumann, prescribed derivative
%   (c) r.boundaries{d}(i,j)  = 1 gives Dirichlet, prescribed field
%   Boundary values are set on field i through boundvalue(i,..,j,..).
%   This array is dynamically returned by r.boundfun{dir}(a,r)
%   and are set equal to interior derivatives to give an approximate value
%   First dimension is the field index, last dimension is the ensemble
%   If r.indext = 0, there is no time index, otherwise there is one.
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License

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
  elseif r.boundaries{dir}(i,1) == -1    %% Use Neuman at low end?
    d2(i1,:,1,:) = 2*((a(i,:,2,:)-a(i,:,1,:))/dx-bval(i,:,1,:))/dx;
  else                                   %% Derivative is periodic
    d2(i1,:,1,:) = (a(i,:,2,:)-2*a(i,:,1,:)+a(i,:,e,:))/dx2;
  end                                    %% End lower boundary case
  if r.boundaries{dir}(i,2)  == 1        %% Use Dirichlet at high end?
    d2(i1,:,e,:) = (a(i,:,e-1,:)-bval(i,:,2,:))/dx2;
  elseif  r.boundaries{dir}(i,2)  == -1  %% Use Neuman at high end?
    d2(i1,:,e,:) = 2*(bval(i,:,2,:)-(a(i,:,e,:)-a(i,:,e-1,:))/dx)/dx;
  else                                   %% Derivative is periodic
    d2(i1,:,e,:) = (a(i,:,1,:)-2*a(i,:,e,:)+a(i,:,e-1,:))/dx2;
  end                                    %% End upper boundary case
i1 = i1+1;                               %% Increment derivative index
end                                      %% End loop over field indices
d2 = reshape(d2,sz);                     %% reshape to output matrix size
end                                      %% end xd2 function