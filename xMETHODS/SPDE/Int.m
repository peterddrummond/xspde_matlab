function o  =  Int(o,varargin)                 
%   a  =  INT(o,[dx,c, bounds,] p) integrates an input over a lattice. 
%   Input:  field 'o', last input are the parameters in structure 'p'.
%   Optional inputs: stepsize 'dx', cell 'c', integration limits 'bounds'.
%   Note: 'dx' should have total dimension equal to p.dimensions;
%   c= 1 is the default if none is input: it gives the default field shape;
%   'bounds' is an array of size [p.dimensions,2],  which specifies 
%   the lower and upper integration bounds in each direction.
%   Output: space integral of 'o' returned in reduced lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   If p.indext = 1 can integrate over time, otherwise not
%   Dimensions above p.dimensions are not integrated.
%   If no dx is present, p.dx is assumed, integrates over all space.
%   First dimension is the field index or line index
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

c = 1;                                           %% Default cell
switch nargin                                    %% Check for arguments
    case 5                                       %% Case of four arguments
        dx = varargin{1};                        %% dx is input
        c =  varargin{2};                        %% cell index is input
        b =  varargin{3};                        %% limits are input
        p =  varargin{4};                        %% p is the last argument
    case 4                                       %% Case of four arguments
        dx = varargin{1};                        %% dx is input
        c =  varargin{2};                        %% cell index is input
        p =  varargin{3};                        %% p is the last argument
    case 3                                       %% Case of three arguments
        dx =varargin{1};                         %% dx is input
        p = varargin{2};                         %% p is the last argument
        b = zeros(p.dimensions,2);               %% bounds are set to zeros
    case 2                                       %% Case of two arguments
        p = varargin{1};                         %% p is second argument
        dx = p.dx;                               %% dx is the step 
        dx(1) = 0;                               %% remove time step
        b = zeros(p.dimensions,2);               %% bounds set to zeros
    otherwise                                    %% Throw an error message
        error ('xSPDE error: Int requires 2, 3 or 4 input arguments.')
end                                              %% End checks on arguments
o = xshape(o,c,p.indext,p);                      %% Reshape to standard
maxind = numel(size(o))- p.indext;               %% Get input dimension
maxind = min(maxind,length(dx));                 %% Get maximum dimension
dx(1) = p.indext*dx(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INTEGRATE SPACE DIMENSIONS 

for i = 1:maxind                                 %% Loop over dimension
    so = size(o);                                %% Get current size
    d = dx(i);                                   %% Get step-size 
    x = p.xc{i};                                 %% Get coordinate vector
    ind = i+p.indext;                            %% Integration index 
    if d > 0                                     %% If integration needed 
      if i>numel(b(:,1)) || b(i,1)>=b(i,2)       %% if no limits
         o = sum(o,ind)*d;                       %% Integrate over all x(i)
      else                                       %% Integrate with limits
         i1 = prod(so(1:(ind-1)));               %% Get indices below
         i2 = prod(so((ind+1):end));             %% Get indices above 
         o = reshape(o,[i1,so(ind),i2]);         %% Reshape for integration
         o1 = zeros([i1,1,i2]);                  %% Initialize integral
         so(ind) = 1;                            %% Integrated index --> 1 
         j= 1 + floor((b(i,1)-x(1))/d);          %% Get first index
         while x(j) <= b(i,2)                    %% While index is in bound
           o1(:,1,:) = o1(:,1,:)+o(:,j,:)*d;     %% Integrate o over dx
           j = j+1;                              %% Check last index
         end                                     %% End while in bound
         o = reshape(o1,so);                     %% Reshape 
      end                                        %% End if invalid limits
    end                                          %% End if i-th integration
end                                              %% End loop on dimension 
end                                              %% End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XINT FUNCTION