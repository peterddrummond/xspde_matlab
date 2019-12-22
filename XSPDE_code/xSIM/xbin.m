function bin  =  xbin(o,varargin)                 
%   a  =  XBIN(o,[dx,] r) bins an input observable over a lattice. 
%   Input:  one component  variable 'o', bin size `dx', lattice 'r'.
%   Note 'dx' should have dimension 1 or at most equal to r.dimension
%   Output: '1' for coordinate inside bin, '0' otherwise.
%   Bins in first direction where dx(i) > 0 or uses dx as a scalar.
%   If no dx is present, in.dx(2)*in.steps(2) is assumed, ie, x-bins.
%   Use this to compute probabilities of average positions
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

switch nargin                                %% checks for arguments
    case 3                                   %% case of three arguments
        dx =varargin{1};                     %% dx is input
        r = varargin{2};                     %% r is the last argument
    case 2                                   %% case of two arguments
        r = varargin{1};                     %% r is the second argument
        dx = r.dx.*r.steps;                  %% dx is the volume element
    otherwise                                %% throw an error message
        error ('xSPDE error: xbin requires 2 or 3 input arguments.')
end                                          %% end checks for arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RESHAPE DATA FOR BINNING  

ldx = length(dx);                            %% Dimension of volume element
sz = size(o);                                %% Get original size to return
if (length(sz) == 2)                         %% If standard matrix input
        o = reshape(o,[sz(1),r.d.int]);      %% Unflatten lattice
end                                          %% End if standard matrix
maxind = min(ndims(o)-1,ldx);                %% Get maximum dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET DIMENSION FOR BINNING

nd  = 2;                                     %% Initial dimension of bin
if ldx >1                                    %% If dx is a vector
    while (dx(nd) <= 0 && nd < maxind)       %% Loop over dx dimension
        nd=nd+1;                             %% Increment dx dimension
    end                                      %% End loop over dx dimension
    dx= dx(nd);                              %% Set dx as a scalar
end                                          %% End if dx is a vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BIN THE DATA OVER ND

if  dx <= 0                                  %% throw error on negative bin
	error ('xSPDE error: xbin requires dx>0')
end                                          %% End error on negative bin
bin1  = r.r{nd}-dx/2;                        %% bin 1 is lower bin boundary
bin2  = r.r{nd}+dx/2;                        %% bin 2 is upper bin boundary
bin   = (o>bin1&o<bin2)/dx;                  %% returns 1/dx if o in bin
end                                          %% end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XBIN FUNCTION