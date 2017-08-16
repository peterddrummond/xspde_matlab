function d  =  xd2(o,varargin) 
%   a = XD2(o,dir,r) calculates a second spatial derivative of a field  
%   uses finite differences with boundary conditions given in in.boundary
%   Boundary conditions are defined at either end: (1,2) sequentially
%   (a) in.boundaries{dir}(1)   = 0 gives the default, periodic 
%   (b) in.boundaries{dir}(1,2) = -1 gives Neumann vanishing derivative
%   (c) in.boundaries{dir}(1,2) = 1 gives Dirichlet, vanishing field 
%   Input: field a, [derivative direction,] lattice r.
%   Output: derivative  d. 
%   If no input direction specified, 
%   a = XD2(o,r) takes derivatives in the x direction, dir=2
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
if nargin < 2 || nargin >3
    error('Function xd2 takes two or three arguments, not %d', nargin)
end
sz = size(o);                            %%get original size to return
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
shape = [1,1,1];
shape(1) =  sz(1)*prod(r.d.int(1:dir)); 
shape(2) =  r.d.int(dir+1); 
shape(3) =  prod(r.d.int(dir+2:end)); 
o = reshape(o,shape);                    %%Unflatten lattice
d = zeros(shape);
dx2 = r.dx(dir)^2;
d(:,2:end-1,:) = (o(:,3:end,:)+o(:,1:end-2,:)-2.*o(:,2:end-1,:))/dx2;
if r.boundaries(dir,1) == 0                %%test for periodic case
   d(:,1,:)   = (o(:,2,:)+o(:,end,:)-2.*o(:,1,:))/dx2;
   d(:,end,:) = (o(:,1,:)+o(:,end-1,:)-2.*o(:,end,:))/dx2;
else
   d(:,1,:)   = (o(:,2,:)-o(:,1,:))/dx2;
   d(:,end,:) = (o(:,end-1,:)-o(:,end,:))/dx2;
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd2 function