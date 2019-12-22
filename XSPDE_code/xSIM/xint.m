function o  =  xint(o,varargin)                 
%   a  =  XINT(o,[dx,] r) integrates an input observable over a lattice. 
%   Input:  one component  variable 'o', measure `dx', lattice 'r'.
%   Note 'dx' should have total dimension equal to r.dimension
%   Output: space integral of 'o' returned in all lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   Last dimension is for ensembles, and is not integrated.
%   If no dx is present, in.dx is assumed to apply, integrates all space.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

switch nargin                                %% Check for arguments
    case 3                                   %% Case of three arguments
        dx =varargin{1};                     %% dx is input
        r = varargin{2};                     %% r is the last argument
    case 2                                   %% Case of two arguments
        r = varargin{1};                     %% r is the second argument
        dx = r.dx;                           %% dx is the volume element
    otherwise                                %% Throw an error message
        error ('xSPDE error: xint requires 2 or 3 input arguments.')
end                                          %% End checks for arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RESHAPE DATA FOR INTEGRATION  

sz = size(o);                                %% Get original size to return
if (length(sz) == 2)                         %% If standard matrix input
        o = reshape(o,[sz(1),r.d.int]);      %% Unflatten lattice
end                                          %% End if standard matrix
maxind = min(ndims(o)-1,length(dx));         %% Get maximum dimension
sizeo = size(o);                             %% Get reshaped size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INTEGRATE SPACE DIMENSIONS 

for i = 2:maxind                             %% Loop over space dimension 
    if dx(i)>0                               %% If i-th integration needed 
        index = ones(1,ndims(o));            %% Initialize index vector
        index(i+1) = sizeo(i+1);             %% Set index to dimension size
        o1 = sum(o,i+1)*dx(i);               %% Integrate o over dx 
        o =repmat (o1,index);                %% Repeat mean over lattice 
    end                                      %% End if i-th integration
end                                          %% End loop to max dimension 
o  = reshape(o,sz);                          %% Restore original dimension
end                                          %% End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XINT FUNCTION