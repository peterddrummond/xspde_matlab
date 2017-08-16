function d  =  xd1(o,varargin) 
%   a = XD1(o,dir,r) calculates a first spatial derivative of a field component 
%   Boundary conditions are defined at either end: (1,2) sequentially
%   (a) in.boundaries{dir}(1)   = 0 gives the default, periodic 
%   (b) in.boundaries{dir}(1,2) = -1 gives Neumann vanishing derivative
%   (c) in.boundaries{dir}(1,2) = 1 gives Dirichlet, vanishing field 
%   Input: field a, [derivative direction,] lattice r.
%   Output: derivative  d. 
%   If no input direction specified, takes derivative in x direction
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
sz = size(o) ;                           %%get original size to return
if nargin == 2
    r = varargin{1};                     %% r is the second argument
    dir = 2;                             %% takes x-derivative
else
    dir =varargin{1};                    %% derivative direction is input
    r = varargin{2};                     %% r is the last argument
end
shape = [1,1,1];
shape(1) =  sz(1)*prod(r.d.int(1:dir)); 
shape(2) =  r.d.int(dir+1); 
shape(3) =  prod(r.d.int(dir+2:end)); 
o = reshape(o,shape);                    %%Unflatten lattice
d = zeros(shape);
dx1 = 2*r.dx(dir);
d(:,2:end-1,:) = (o(:,3:end,:)-o(:,1:end-2,:))/dx1;
if r.boundaries(dir,1) == 0                %%test for periodic case
   d(:,1,:) = (o(:,2,:)-o(:,end,:))/dx1;
   d(:,end,:) = (o(:,1,:)-o(:,end-1,:))/dx1;
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd1 function