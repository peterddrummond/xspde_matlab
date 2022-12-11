function d2  =  xd2(a,varargin) 
%   d2 = XD2(o,dir,r) calculates a second spatial derivative of any field a  
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
                                                            
if nargin < 2 || nargin > 3              %% check number of arguments input
    error('Function xd2 takes two or three arguments, not %d', nargin)
end                                      %% end check number of arguments
if nargin == 2                           %% check if one extra argument
    r = varargin{1};                     %% r is always the last argument
    dir = 2;                             %% x-derivative is default value
else                                     %% two extra arguments input
    dir =varargin{1};                    %% direction is second argument
    r = varargin{2};                     %% r is the last argument
end                                      %% end check if one extra argument
in = r.indext;                           %% get time-index switch
if dir < 2-in || dir > r.dimension       %% check direction input
    error('Derivative direction %d in xd2 is not supported',dir)
end                                      %% end check direction input
sz = size(a);                            %% get input size to return
sz2 = prod(sz(2+in:end));                %% calculate space lattice size
if r.nlattice ~= sz2                     %% check input lattice size
    error('xd2 expects lattice size of %d, not %d',r.nlattice,sz2)
end                                      %% end check input lattice size
bv = r.boundfun(a,dir,r);                %% Call boundary value function
s(1) = sz(1);                            %% Field index
s(2) = prod(sz(2:in+dir-1));             %% Low  indices
s(3) = sz(in+dir);                       %% Derivative index 
s(4) = prod(sz(in+dir+1:end));           %% High indices 
a = reshape(a,s);                        %% Unflatten lattice
d2 = zeros(s);                           %% Initialize derivative to zero
dx = r.dx(dir);                          %% Set the space step value
dx2 = dx^2;                              %% Set the space step squared
e = s(3);                                %% Set the extent
s(3) =  2;                               %% Set shape for boundary values
bv = reshape(bv,s);                      %% Reshape the boundary value
d2(:,:,2:end-1,:) = (a(:,:,3:e,:)+a(:,:,1:e-2,:)-2.*a(:,:,2:e-1,:))/dx2;
for i=1:sz(1)                            %% Loop over field indices
if r.boundaries{dir}(i,1) == 1           %% Use Dirichlet at low end?
    d2(i,:,1,:) = (a(i,:,2,:)-bv(i,:,1,:))/dx2;
elseif r.boundaries{dir}(i,1) == -1      %% Use Neuman at low end?
    d2(i,:,1,:) = 2*((a(i,:,2,:)-a(i,:,1,:))/dx-bv(i,:,1,:))/dx;
else                                     %% Derivative is periodic
    d2(i,:,1,:) = (a(i,:,2,:)-2*a(i,:,1,:)+a(i,:,e,:))/dx2;
end                                      %% End lower boundary case
if r.boundaries{dir}(i,2)  == 1          %% Use Dirichlet at high end?
    d2(i,:,e,:) = (a(i,:,e-1,:)-bv(i,:,2,:))/dx2;
elseif  r.boundaries{dir}(i,2)  == -1    %% Use Neuman at high end?
    d2(i,:,e,:) = 2*(bv(i,:,2,:)-(a(i,:,e,:)-a(i,:,e-1,:))/dx)/dx;
else                                     %% Derivative is periodic
    d2(i,:,e,:) = (a(i,:,1,:)-2*a(i,:,e,:)+a(i,:,e-1,:))/dx2;
end                                      %% End upper boundary case
end                                      %% End loop over field indices
d2 = reshape(d2,sz);                       %% reshape to input matrix size
end                                      %% end xd2 function