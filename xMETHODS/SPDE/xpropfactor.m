function [propagator,da] = xpropfactor (c,p)                       
%   XPROPFACTOR(p)  returns the interaction picture propagation factor.
%   Input:   structure 'p' with momentum cells p.k.
%   Output: 'propagator' is a k-space or constant factor used by xprop.
%           'da' is the parameter needed in N-N boundary corrections
%
%   The propagator is calculated spectrally using the "linear" function.
%   This uses derivatives as p.Dx.. OR p.D{1}.., calculated from p.k{c,d}.
%   The first propagator dimension is either the field index or 1
%   With non-periodic boundaries, ONLY even order derivatives are allowed.
%
%   Note: 'da' is required for Neumann-Neumann spectral boundaries in xprop
%          The value of 'da' for dimension 'd' is obtained from "linear". 
%          ONLY treats 2nd order derivatives in 'd' for the indexed term.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Needs:      p.dtr, p.ipsteps, p.dimensions, p.d.space, p.ranges, p.k
%   Calls:      p.linear{c}(p)
%   Called by:  xensemble
%   Licensed by P. D. Drummond (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  p.D = cell(1,max(4,p.dimensions));            % make cell of derivatives
  dt = p.dtr/p.ipsteps;                         % get effective time step 
  da = 0;                                       % initialise da output
  if p.dimensions > 1                           % if pde case
    for d = 2:p.dimensions                      % loop over dimensions
       p.D{d} = 1i*p.k{c,d};                    % make spectral derivative
    end                                         % end loop on dimensions
    [p.Dx,p.Dy,p.Dz] = deal(p.D{2:4});          % make  Dx,Dy,Dz grids
  end                                           % end if pde
  L =  p.linear{c}(p);                          % get the linear term
  if isequal(L, 0) || isequal(L, [])            % if L zero or null
    propagator = 1;                             % set propagator to 1
  else                                          % if L not zero or null
    sz = size(L);                               % get size of L
    if p.dimensions > 1                         % if pde case
      L = reshape(L,[sz(1),p.d.space{c}]);      % reshape to full size
      if p.setboundaries                        % if setting boundaries?
        da = zeros(p.fields{c},p.dimensions);   % Neumann-Neumann da term
        fda = dt/pi^2;                          % compute factor for da
        for d = 2:p.dimensions                  % loop over dimensions
          in = num2cell(ones(1,p.dimensions-1));% get cell of indices = 1
          L1 = L(:,in{:});                      % get L of 1st index in 'd'
          in{d-1} = 2;                          % get 2nd index in 'd'
          L2 = L(:,in{:});                      % get L of 2nd index in 'd'
          da(:,d) = -(L2-L1)*fda*p.ranges(d)^2; % compute NN correction
        end                                     % end loop over dimensions
      end                                       % end if setting boundaries
    end                                         % end if pde
    propagator = exp(L*dt);                     % get factor
  end                                           % end if L nonzero
end                                             % end xpropfactor