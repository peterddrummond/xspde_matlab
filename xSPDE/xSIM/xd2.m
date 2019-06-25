function d  =  xd2(a,varargin) 
%   a = XD2(o,dir,r) calculates a second spatial derivative of any field a  
%   using second order (in dx) central finite difference methods.  
%   Output: second derivative  d.
%   Boundary conditions are in the field, except for periodic case
%   Derivatives at boundaries are for measurements, not for propagation
%   and are set equal to interior derivatives to give an approximate value
%   If two inputs, takes derivative in x direction
%   If three inputs, takes transverse derivative in axis  'dir' = 2,3..
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License 
                                                            
if nargin < 2 || nargin >3
    error('Function xd2 takes two or three arguments, not %d', nargin)
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
    error('Derivative direction %d in xd2 is not supported',dir)
end
sz2 = prod(sz(2:end));
if r.nlattice ~= sz2
    error('xd2 expects lattice size of %d, not %d',r.nlattice,sz2)
end
ds=dir-1;
shape(1) =  sz(1);                       %%Field index
shape(2) =  prod(r.d.int(1:ds));      %%Low space indices
shape(3) =  r.d.int(dir);                %%Space derivative index 
shape(4) =  prod(r.d.int(dir+1:end));    %%High space/ensemble indices 
a = reshape(a,shape);                    %%Unflatten lattice
d = zeros(shape);
dx = r.dx(dir);
dx2 = dx^2;
shape(3) =  2;                            %%for boundary values
cellboundvalue=r.boundfun(a,r);
boundvalue = reshape(cellboundvalue{dir},shape);
d(:,:,2:end-1,:) = (a(:,:,3:end,:)+a(:,:,1:end-2,:)-2.*a(:,:,2:end-1,:))/dx2;
for i=1:sz(1)
if r.boundaries{dir}(i,1) == 1
    d(i,:,1,:) = d(i,:,2,:);
elseif r.boundaries{dir}(i,1) == -1
    d(i,:,1,:) = 2*((a(i,:,2,:)-a(i,:,1,:))/dx2-boundvalue(i,:,1,:)/dx);
else
    d(i,:,1,:) = (a(i,:,2,:)-2*a(i,:,1,:)+a(i,:,end,:))/dx2;       %%periodic case
end
if r.boundaries{dir}(i,2)  == 1
    d(:,:,end,:) = d(:,:,end-1,:);
elseif  r.boundaries{dir}(i,2)  == -1
    d(:,:,end,:) = 2*(boundvalue(i,:,end,:)/dx-(a(:,:,end,:)-a(:,:,end-1,:))/dx2);
else
    d(:,:,end,:) = (a(:,:,1,:)-2*a(:,:,end,:)+a(:,:,end-1,:))/dx2; %%periodic case
end
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd2 function