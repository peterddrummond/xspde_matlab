function o  =  xint(o,varargin)                 
%   a  =  XINT(o,[dx,] r) integrates an input observable over a lattice. 
%   Input:  one component  variable 'o', measure `dx', lattice 'r'.
%   Note 'dx' should have total dimension equal to r.dimension
%   Output: space integral of 'o' returned in all lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   Second dimension is for ensembles, and is not integrated.
%   If no dx is present, in.dx is assumed to apply, integrates all space.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

switch nargin                                %% checks for arguments
    case 3                                   %% case of three arguments
        dx =varargin{1};                     %% dx is input
        r = varargin{2};                     %% r is the last argument
    case 2                                   %% case of two arguments
        r = varargin{1};                     %% r is the second argument
        dx = r.dx;                           %% dx is the volume element
    otherwise                                %% throw an error message
        error ('xSPDE error: xint requires 2 or 3 input arguments.')
end                                          %% end checks for arguments
sz = size(o);                                %%get original size to return
if (length(sz) == 2)                         %%Standard matrix input
        o = reshape(o,[sz(1),r.d.int]);      %%Unflatten lattice
end
maxind = min(length(o)-2,length(dx));        %% get maximum dimension
sizeo = size(o);
for i = 2:maxind                             %% loop over space dimension 
    if dx(i)>0                               %% is i-th integration needed? 
        index = ones(1,length(o));           %% Initialize index vector
        index(i+2) = sizeo(i+2);             %% Set index to dimension size
        %o1 = trapz(o,i+2);
        o1 = sum(o,i+2);
        o =repmat (o1*dx(i),index);          %% Repeat mean over lattice 
    end                                      %% end i-th integration
end                                          %% end loop to max dimension 
o  = reshape(o,sz);                          %restore original dimension
end                                          %%end function