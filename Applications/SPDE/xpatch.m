function [a,pa,sw,lp]  =  xpatch(ac,c,p) 
%   [a,pa,sw,lp] = XPATCH(ac,c,p) calculates 1D spatial patches
%   for trigonometric spectral transform methods used in cell c.
%   Returns  corrected field a, the patch, switch sw, transformed patlp
%   Here patlp is the patch Laplacian multiplied by dt*D_j.
%   If sw = 0, use patlp; if sw = 1, calculate numerically
%   The field is assumed periodic unless there are specified boundaries
%   A patch is returned in any dimension, but with only one nonzero patch
%   Boundary types are defined through p.boundaries{c,d}(i,b) for cell c,
%   space dimension (d=2,3..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d>1 is the transverse space dimension index, and options are:
%
%   (a) p.boundaries{c,d}(i,b)  = 0 gives the default, periodic 
%   (b) p.boundaries{c,d}(i,b)  = -1 gives Neumann, prescribed derivative
%   (c) p.boundaries{c,d}(i,b)  = 1 gives Dirichlet, prescribed field
%
%   Boundary values are specified through p.boundfun{a,c,d,p}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: xprop,MPx
%   Needs: p.setboundaries, p.propagator, p.dimensions, p.boundaries
%   Licensed by P. D. Drummond (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sw = 0;                                          % initial switch
a = ac{c};                                       % get field in cell c
pa = 0;                                          % initial patch
lp = 0;                                          % initial laplacian term
sx = size(a);                                    % get size
irange = 1:sx(1);                                % get range of field index
for d = 2:p.dimensions                           % loop on space dimensions
  bnd = p.boundaries{c,d};                       % store boundary type
  type = bnd(irange,1)+2*bnd(irange,2);          % boundary type index    
  if any(type) && sx(d) > 1                      % If index not periodic
    i1 = p.ind;                                  % initialize index cells
    pa = zeros(sx);                              % initial patch
    lp = zeros(sx);                              % initial laplacian term
    b  = p.boundfun(ac,c,d,p);                   % get boundary values                      
    xa = p.origins(d);                           % get origin
    dx = p.ranges(d);                            % get range
    x  = p.r{d} - xa;                            % relative x coordinate
    for i = irange                               % loop over components
      i1{1} = i;                                 % set first index to i
      switch type(i)                             % switch on boundary case
        case -3                                  % Robin-Robin boundaries
          db = (b{i,2} - b{i,1})/dx;             % difference in boundaries
          pa(i1{:}) = pa(i1{:}) + b{i,1}.*x + 0.5*x.^2.*db;  % patch
          lp(i1{:}) = lp(i1{:})+db*p.da{c,d}(i); % dynamic bdary correction
        case -1                                  % Dirichlet-Robin boundary
          pa(i1{:}) = pa(i1{:}) + b{i,1} + x.*b{i,2}; % patch
        case 0                                   % periodic in [c,d,i]!
        case 1                                   % Robin-Dirichlet
          pa(i1{:}) = pa(i1{:}) + b{i,2} + (x - dx).*b{i,1}; % patch
        case 3                                   % Dirichlet-Dirichlet
          db = (b{i,2} - b{i,1})/dx;             % difference in boundaries          
          pa(i1{:}) = pa(i1{:}) + b{i,1} + x.*db;% Dirichlet-Dirichlet
        otherwise                                % invalid boundary type
          error('Invalid boundary: [%d,%d]',bnd(i,1),bnd(i,2));
      end                                        % end switch     
    end                                          % end irange loop
  end                                            % end if type not periodic
end                                              % end dimension loop
a = a - pa;                                      % return corrected field
end                                              % end xpatch function                                        