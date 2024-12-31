function Da  =  DS(a,varargin) 
%   d1 = DS(a,[n,d,c,j],p) calculates n-th spatial derivatives
%   using spectral transform methods, either Fourier or trigonometric  
%   Output: the n-th derivative in dimension d, cell c, indices j. Returns:
%
%   (i)   2 inputs: 2nd derivatives,  x-dimension, cell 1, all indices  
%   (ii)  3 inputs: n-th derivatives, x-dimension, cell 1, all indices  
%   (ii)  4 inputs: n-th derivatives, dimension d, cell 1, all indices
%   (iii) 5 inputs: n-th derivatives, dimension d, cell c, all indices 
%   (iii) 6 inputs: n-th derivatives, dimension d, cell c, index list j
%   
%   If j is input as a vector, indices in the list are returned in order.
%   In this case, the output derivative may have a different size to a.
%   The field is assumed periodic unless there are specified boundaries
%   Boundary types are defined through p.boundaries{c,d}(i,b) for cell c,
%   space dimension (d=2,3..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d=2,3.. is the transverse space dimension (x,y,..). Options are:
%
%   (a) p.boundaries{c,d}(i,b)  = 0 gives the default, periodic 
%   (b) p.boundaries{c,d}(i,b)  = -1 gives Neumann, prescribed derivative
%   (c) p.boundaries{c,d}(i,b)  = 1 gives Dirichlet, prescribed field
%
%   Boundary values are specified through p.boundval{c,d}{i,b}
%   This requires p.boundval to be updated before it is used.
%   Boundary values match the field except for having only one first index
%   In the dimension of the derivative (d), one also has a single index 
%   Note that b = 1,2 specify values at the upper and lower boundaries  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: deriv, observe, output
%   Needs: p.setboundaries, p.propagator, p.dimensions, p.boundaries
%   Calls: xdct1,xdct2,xdct3,xdst1,xdst2,xdst3,fft,ifft
%   Licensed by P. D. Drummond, S. Kiesewetter (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sx = size(a);                                    % get original size
irange = 1:sx(1);                                % get range of all indices
p = varargin{end};                               % get p: the last argument
n = varargin{1};                                 % get order if input
d = 2;                                           % set default dimension
c = 1;                                           % set default cell
switch nargin
  case 2
    n = 2;                                       % second derivative
  case 3
  case 4
	d = varargin{2};                             % component list is input 
  case 5 
    d = varargin{2};                     
	c = varargin{3};                             % component list is input
  case 6 
    d = varargin{2};   
	irange = varargin{4};                        % component list is input
    c = varargin{3};                             % cell index is input
  otherwise
    error('Function DS takes two to six arguments, not %d', nargin)
end
irl  = length(irange);
sz   = size(a);                                  % field size
s1   = sz(1);                                    % internal index size
b    = num2cell(zeros(s1,2));                    % initial boundary value
bnd  = p.boundaries{c,d};                        % store boundary type
type = bnd(irange,1)+2*bnd(irange,2);            % boundary type index
if any(type)                                     % If an index not periodic
    b  = p.boundval{c,d};                        % get boundary values                      
    xa = p.origins(d); dx = p.ranges(d);         % get origin and range
    x  = p.r{d} - xa;                            % relative x coordinate
    switch n
      case 1
        DC = [1,x];
      case 2
        DC = [0,1];
      otherwise
        DC = [0,0];
    end  
end                                              % end if not periodic
j = 1;
szD = [irl,sz(2:end)];
Da = zeros(szD);
for i = irange                                   % loop over components
    in = p.ind;                                  % get input indices
    out = in;                                    % get output indices
    in{1} = i;                                   % input field index = i
    out{1} = j;                                  % output field index = i
    et = 0;                                      % initial dynamic bdary
    switch type(j)                               % switch on boundary case
        case -3                                  % Robin-Robin boundaries
          db = (b{i,2} - b{i,1})/dx;             % difference in boundaries
          b1 = b{i,1}.*x + 0.5*x.^2.*db;         % static bdary correction
          et = DC(1)*b{i,1}+DC(2)*db;            % dynamic bdary correction
          D  = a(in{:}) - b1;                    % subtract correction
          D  = xdct1(D,d);                       % cosine-1 transform
        case -1                                  % Dirichlet-Robin boundary
          b1 = b{i,1} + x.*b{i,2};               % static bdary correction
          et = DC(1)*b{i,2};                     % dynamic bdary correction
          D  = a(in{:}) - b1;                    % subtract correction           
          D  = xdst3(D,d);                       % take sine-3 transform
        case 0                                   % periodic boundaries
          D  = fft(a(in{:}),[],d);               % Fourier transform
        case 1                                   % Robin-Dirichlet
          b1 = b{i,2} + (x - dx).*b{i,1};        % static bdary correction
          et = DC(1)*b{i,1};                     % dynamic bdary correction
          D  = a(in{:}) - b1;                    % subtract correction
          D  = xdct3(D,d);                       % take cosine-3 transform
        case 3                                   % Dirichlet-Dirichlet
          db = (b{i,2} - b{i,1})/dx;             % difference in boundaries          
          b1 = b{i,1} + x.*db;                   % Dirichlet-Dirichlet
          et = DC(1)*db;                         % dynamic bdary correction
          D  = a(in{:}) - b1;                    % subtract correction
          D  = xdst1(D,d);                       % take sine-1 transform            
        otherwise                                % invalid boundary type
          error('Invalid boundary: [%d,%d]',bnd(i,1),bnd(i,2));
    end                                          % end switch     
    D = ((1i*p.k{c,d}(in{:})).^n).*D;            % differentiate in k-space
    switch type(j)*(-1)^n                        % inverse type = t*(-1)^n
        case -3                                  % Robin-Robin boundaries
          D = xdct1(D,d);                        % cosine-1 transform
        case -1                                  % Dirichlet-Robin
          D = xdst2(D,d);                        % sine-2 transform
        case 0                                   % Dirichlet-Robin
          D = ifft(D,[],d);                      % Fourier transform
        case 1                                   % Robin-Dirichlet 
          D = xdct2(D,d);                        % cosine-2 transform
        case 3                                   % Dirichlet-Dirichlet       
          D = xdst1(D,d);                        % sine-1 transform
    end                                          % end switch
    Da(out{:}) = D + et;                         % store derivative
    j = j + 1;                                   % increment j index
end                                              % end components loop
end                                              % end DS function                                          