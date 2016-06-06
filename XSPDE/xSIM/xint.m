function o  =  xint(o,varargin)                 
%   a  =  XINT(o,[dx,] r) integrates an input observable over a lattice. 
%   Input:  one component  variable 'o', measure `dx', lattice 'r'.
%   Note 'o' and 'dx' should have total dimension equal to r.dimension
%   Output: space integral of 'o' returned in all lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   First observable dimension is for ensembles, and is not integrated.
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
maxind = min(r.dimension,length(dx));        %% get maximum dimension
inputsize = size(o);                         %% input size 
o = reshape(o,r.d.int);                      %% unflatten lattice 
%xintsize = size(o)
for i = 2:maxind                             %% loop over space dimension 
    if dx(i)>0                               %% is i-th integration needed? 
        index = ones(1,r.dimension+1);       %% set index to ones 
        index(i+1) = r.d.int(i+1);           %% set i-th index to size
        o =repmat(sum(o,i+1)*dx(i),index);   %% integrate and repeat
    end                                      %% end i-th integration
end                                          %% end loop to max dimension 
o  = reshape(o,inputsize);                   %Flatten lattice 
end                                          %%end function