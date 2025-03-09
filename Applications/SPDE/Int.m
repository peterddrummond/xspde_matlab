function o  =  Int(o,varargin)                 
%   a  =  INT(o,[dx, bounds, c] p) integrates an input over a lattice. 
%   Input:  field 'o', last input are the parameters in structure 'p'.
%   Optional inputs: stepsize 'dx', integration limits 'bounds', cell 'c'.
%   Note: 'dx' should have total dimension equal to p.dimensions;
%   c= 1 is the default if none is input: it gives the default field shape;
%   This is only required if the input field needs reshaping
%   'bounds' is an array of size [p.dimensions,2],  which specifies 
%   the lower and upper integration bounds in each direction.
%   If omitted, all space points are integrated
%   Output: space integral of 'o' returned in reduced lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   If p.indext = 1, can integrate over time, (in output functions only)
%   Otherwise there is no time integration.
%   Dimensions above p.dimensions are not integrated.
%   If no dx is present, p.dx is assumed, integrates over all space.
%   Note: if transformed fields are input, dk must be included.
%   First dimension is the field index or line index
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

p  = varargin{end};                              %% p is the last argument
dx = p.dx;                                       %% Default dx 
b  = zeros(p.dimensions,2);                      %% Default bounds
c  = 1;                                          %% Default cell
switch nargin                                    %% Check for arguments
    case 5                                       %% Case of four arguments
        dx = varargin{1};                        %% dx is input
        b =  varargin{2};                        %% limits is input
        c =  varargin{3};                        %% cell index are input
    case 4                                       %% Case of four arguments
        dx = varargin{1};                        %% dx is input
        b =  varargin{2};                        %% limits are input
    case 3                                       %% Case of three arguments
        dx = varargin{1};                        %% dx is input
    case 2                                       %% Case of two arguments
    otherwise                                    %% Throw an error message
      error ('xSPDE error: Int requires 2, 3, 4 or 5 input arguments.')
end                                              %% End checks on arguments
o = xshape(o,c,p.indext,p);                      %% Reshape to standard
maxind = numel(size(o))- p.indext;               %% Get input dimension
maxind = min(maxind,length(dx));                 %% Get maximum dimension
dx(1)  = p.indext*dx(1);

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