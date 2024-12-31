function [propagator,da,propx] = xpropfactor(p)                       
%   XPROPFACTOR(p)  returns the interaction picture propagation factors.
%   Input:   structure 'p' with momentum cells p.k.
%   Output: 'propagator{c}' is a k-space or constant factor used by xprop
%           'da{c,d}' is the parameter needed in N-N boundary corrections
%           'propx{c}' is a second k-space factor used by MPx
%
%   The propagator is calculated spectrally using the "linear" function.
%   This uses derivatives as p.Dx.. OR p.D{1}.., calculated from p.k{c,d}.
%   The first propagator dimension is either the field index or 1
%   With non-periodic boundaries, ONLY even order derivatives are allowed.
%
%   Note: 'da' is required for Neumann-Neumann spectral boundaries in xprop
%          The value of 'da' for dimension 'd' is obtained from "linear". 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Needs:      p.dtr, p.ipsteps, p.dimensions, p.d.space, p.ranges, p.k
%   Calls:      p.linear{c}(p)
%   Called by:  xensemble
%   Licensed by P. D. Drummond (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

da = cell(1,p.fieldcells);                       % initialize da cells
propagator = da;                                 % initialize propagator
propx = da;                                      % initialize propagator
for c = 1:p.fieldcells                           % loop over field cells
  propagator{c} = 1;                             % get default propagator
  dt = p.dtr/p.ipsteps;                          % get effective time step 
  propx{c}  = dt;                                % get default propx
  p.D = cell(1,max(4,p.dimensions));             % make cell of derivatives
  if p.dimensions > 1                            % if pde case
    for d = 2:p.dimensions                       % loop over dimensions
       p.D{d} = 1i*p.k{c,d};                     % make spectral derivative
    end                                          % end loop on dimensions
    [p.Dx,p.Dy,p.Dz] = deal(p.D{2:4});           % make  Dx,Dy,Dz grids
  end                                            % end if pde
  if ~isequal(p.linear,[]) && ~isequal(p.linear{c},0) 
    L = p.linear{c}(p);                          % get the linear term
    propagator{c} = exp(L*dt);                   % get propagator
    for d = 2:p.dimensions                       % loop over dimensions
      da{c,d} = 0;                               % initialise da output
      szL = size(L,1);                           % get size of L
      L = reshape(L,[szL,p.d.space{c}]);         % reshape to full size
      if p.dimensions > 1 && p.setboundaries     % if setting boundaries?
        fda = dt/pi^2;                           % compute factor for da
        in = num2cell(ones(1,p.dimensions-1));   % get cell of indices = 1
        L1 = L(:,in{:});                         % get L of 1st index in d
        in{d-1} = 2;                             % get 2nd index in d
        L2 = L(:,in{:});                         % get L of 2nd index in d
        da{c,d} = -(L2-L1)*fda*p.ranges(d)^2;    % compute NN correction
      end                                        % end if set boundaries
    end                                          % end dimensions loop
    if any(L == 0)
        L = L+1.e-10; 
    end
    propx{c} =  (exp(L*dt)-1)./L;                % get propagator
  end                                            % end if ~isequal
end                                              % end cell loop
end                                              % end xpropfactor