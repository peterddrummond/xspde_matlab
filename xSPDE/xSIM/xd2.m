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
if dir < 2 || dir > r.dimension          %% check direction input
    error('Derivative direction %d in xd2 is not supported',dir)
end                                      %% end check direction input
sz = size(a);                            %% get input size to return
sz2 = prod(sz(2:end));                   %% calculate input lattice size
if r.nlattice ~= sz2                     %% check input lattice size
    error('xd2 expects lattice size of %d, not %d',r.nlattice,sz2)
end                                      %% end check input lattice size
shape(1) =  sz(1);                       %% Field index
shape(2) =  prod(r.d.int(1:dir-1));      %% Low space indices
shape(3) =  r.d.int(dir);                %% Space derivative index 
shape(4) =  prod(r.d.int(dir+1:end));    %% High space/ensemble indices 
a = reshape(a,shape);                    %% Unflatten lattice
d = zeros(shape);                        %% Initialize derivative to zero
dx = r.dx(dir);                          %% Set the space step value
dx2 = dx^2;                              %% Set the space step squared
shape(3) =  2;                           %% Set shape for boundary values
boundvalue=r.boundfun(a,dir,r);          %% Call boundary value function
boundvalue = reshape(boundvalue,shape);  %% Reshape the boundary value
d(:,:,2:end-1,:) = (a(:,:,3:end,:)+a(:,:,1:end-2,:)-2.*a(:,:,2:end-1,:))/dx2;
for i=1:sz(1)                            %% Loop over field indices
if r.boundaries{dir}(i,1) == 1           %% Use Dirichlet at low end?
    d(i,:,1,:) = d(i,:,2,:);             %% Derivative at low end defined
elseif r.boundaries{dir}(i,1) == -1      %% Use Neuman at low end?
    d(i,:,1,:) = 2*((a(i,:,2,:)-a(i,:,1,:))/dx2-boundvalue(i,:,1,:)/dx);
else                                     %% Derivative is periodic
    d(i,:,1,:) = (a(i,:,2,:)-2*a(i,:,1,:)+a(i,:,end,:))/dx2;
end                                      %% End lower boundary case
if r.boundaries{dir}(i,2)  == 1          %% Use Dirichlet at high end?
    d(i,:,end,:) = d(i,:,end-1,:);       %% Derivative at high end defined
elseif  r.boundaries{dir}(i,2)  == -1    %% Use Neuman at high end?
    d(i,:,end,:) = 2*(boundvalue(i,:,end,:)/dx-(a(i,:,end,:)-a(i,:,end-1,:))/dx2);
else                                     %% Derivative is periodic
    d(i,:,end,:) = (a(i,:,1,:)-2*a(i,:,end,:)+a(i,:,end-1,:))/dx2;
end                                      %% End upper boundary case
end                                      %% End loop over field indices
d = reshape(d,sz);                       %% reshape to input matrix size
end                                      %% end xd2 function