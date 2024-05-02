function [a,varargout]  =  MP(a,xi,p)            
%   a = MP(a,xi,p) propagates a step using the midpoint + IP method.
%   Input:  cell of fields 'a',  noises 'xi', parameter structure 'p'.
%   Output: cell of fields 'a',  with derivative and linear propagator. 
%   Includes jumps if specified by p.jump > 0
%   All derivatives are calculated at the midpoint
%   All boundary values are calculated at the first and last point
%   Licensed by Peter D. Drummond, (2024) - see License and manual.

dt  = p.dtr/2.0;                                % store half time-step
[a,p.boundval] =  p.prop(a,p);                  % propagate to midpoint
p.t = p.t + dt;                                 % get midpoint time
a0 = a;                                         % store initial cells
for iter = 1:p.iterations                       % midpoint iteration loop
  for c = 1:p.fieldcells                        % loop over field cells
    a{c} = a0{c} + reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c}); % derivs
    if p.setboundaries                          % if boundaries
       [a,p.boundval] = p.setbound(a,p);        % set boundary values
    end                                         % end if boundaries
   end                                          % end field cells loop
end                                             % end iteration loop

for c = 1:p.fieldcells                          % loop over field cells
  a{c} = 2*a{c}-a0{c};                          % add midpoint derivatives
end                                             % end loop over field cells
p.t = p.t + dt;                                 % get endpoint time
a = p.prop(a,p);                                % propagate to end point 
varargout = {1,2,2,1,1};                        % method parameters
%% St. order = 1, det. order = 2, ipsteps = 2, vect = 1(Y), cell  = 1(Y)
end                                             %  end function call    