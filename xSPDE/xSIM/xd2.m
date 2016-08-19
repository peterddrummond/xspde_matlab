function d  =  xd2(o,varargin) 
%   a = XD2(o,dir,r) calculates a second spatial derivative of a field  
%   uses finite derivatives with boundary conditions given in in.boundary
%   Here, in.boundaries(dir)=0 is the default, periodic case
%   or, in.boundaries(dir)1 is the Neumann, vanishing derivative case
%   and, in.boundaries(dir)=2 is the alternate Dirichlet, vanishing field case
%   Input: field a, [derivative direction,] lattice r.
%   Output: derivative  d. 
%   If no input direction specified, takes derivative in x direction, dir=2
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
sz = size(o);                            %%get original size to return
if nargin == 2
    r = varargin{1};                     %% r is the second argument
    dir = 2;                             %% x-derivative
else
    dir =varargin{1};                    %% direction is input
    r = varargin{2};                     %% r is the last argument
end
shape = [1,1,1];
shape(1) =  sz(1)*prod(r.d.int(1:dir)); 
shape(2) =  r.d.int(dir+1); 
shape(3) =  prod(r.d.int(dir+2:end)); 
o = reshape(o,shape);                    %%Unflatten lattice
d = zeros(shape);
dx2 = r.dx(dir)^2;
d(:,2:end-1,:) = (o(:,3:end,:)+o(:,1:end-2,:)-2.*o(:,2:end-1,:))/dx2;
if r.boundaries(dir) == 0                %%test for periodic case
   d(:,1,:) = (o(:,2,:)+o(:,end,:)-2.*o(:,1,:))/dx2;
   d(:,end,:) = (o(:,1,:)+o(:,end-1,:)-2.*o(:,end,:))/dx2;
elseif r.boundaries(dir) == 1            %%test for Neumann case
   d(:,1,:) = (o(:,2,:)-o(:,1,:))/dx2;
   d(:,end,:) = -(o(:,end,:)-o(:,end-1,:))/dx2;
end
d = reshape(d,sz);                       %%reshape to matrix
end                                      %%end xd2 function