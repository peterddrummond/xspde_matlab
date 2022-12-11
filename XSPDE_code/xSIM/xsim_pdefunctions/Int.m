function o  =  Int(o,varargin)                 
%   a  =  INT(o,[dx,][bounds,] r) integrates an input over a lattice. 
%   Input:  field 'o', lattice 'r'.
%   Optional inputs:   stepsize 'dx',  integration bounds 'bounds'.
%   Note 'dx' should have total dimension equal to r.dimension,
%   'bounds' is an array of size [r.dimension,2],  which specifies 
%   the lower and upper integration bounds in each direction.
%   Output: space integral of 'o' returned in reduced lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   If r.indext = 1 can integrate over time, otherwise not
%   Any dimensions above r.dimension are not integrated.
%   If no dx is present, in.dx is assumed, integrates over all space.
%   First dimension is the field index or line index
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS


switch nargin                                %% Check for arguments
    case 4                                   %% Case of four arguments
        dx =varargin{1};                     %% dx is input
        b =varargin{2};                      %% bounds is input
        r = varargin{3};                     %% r is the last argument
    case 3                                   %% Case of three arguments
        dx =varargin{1};                     %% dx is input
        r = varargin{2};                     %% r is the last argument
        b = zeros(r.dimension,2);            %% bounds is set to zeros
    case 2                                   %% Case of two arguments
        r = varargin{1};                     %% r is the second argument
        dx = r.dx;                           %% dx is the volume element
        dx(1) = 0;                           %% remove time step
        b = zeros(r.dimension,2);            %% bounds set to zeros
    otherwise                                %% Throw an error message
        error ('xSPDE error: xint requires 2, 3 or 4 input arguments.')
end                                          %% End checks for arguments
o = xshape(o,r.indext,r);                    %% Reshape to xspde standards
maxind = numel(size(o))- r.indext;           %% Get input dimension
maxind = min(maxind,length(dx));             %% Get maximum dimension
dx(1) = r.indext*dx(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INTEGRATE SPACE DIMENSIONS 

for i = 1:maxind                             %% Loop over dimension
    so = size(o);                            %% Get current size
    d = dx(i);                               %% Get step-size 
    x = r.xc{i};                             %% Get coordinate vector
    ind = i+r.indext;                        %% Integration index in field
    if d > 0                                 %% If i-th integration needed 
      if i>numel(b(:,1)) || b(i,1)>=b(i,2)   %% if no limits
         o = sum(o,ind)*d;                   %% Integrate o over all x(i)
      else                                   %% Integrate with limits
         i1 = prod(so(1:(ind-1)));           %% Get indices below
         i2 = prod(so((ind+1):end));         %% Get indices above 
         o = reshape(o,[i1,so(ind),i2]);     %% Reshape for integration
         o1 = zeros([i1,1,i2]);              %% Initialize integral
         so(ind) = 1;                        %% Integrated index --> one  
         j= 1 + floor((b(i,1)-x(1))/d);      %% Get first index
         while x(j) <= b(i,2)                %% While index is in bound
           o1(:,1,:) = o1(:,1,:)+o(:,j,:)*d; %% Integrate o over dx
           j = j+1;                          %% Check last index
         end                                 %% End while index is in bound
         o = reshape(o1,so);                 %% Reshape to xSPDE dimensions
      end                                    %% End if invalid limits
    end                                      %% End if i-th integration
end                                          %% End loop to max dimension 
end                                          %% End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XINT FUNCTION