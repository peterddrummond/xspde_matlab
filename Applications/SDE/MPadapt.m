function [a,varargout]  =  MPadapt(a,xi,p)            
%   ac = MP(ac,xi,p) propagates a step using the midpoint + IP method.   
%   Input:  cell of fields 'a',  noises 'xi', parameter structure 'p'.
%   Output: cell of fields 'a',  with derivative and linear propagator. 
%   Licensed by Peter D. Drummond, (2024) - see License and manual.

[a,p.boundval] =  p.prop(a,p);                   % Propagate the fields
dt  = p.dtr/2.0;                                 % Half time-step 
p.t = p.t + dt;                                  % get midpoint time
inv = cell(1,p.fieldcells); ac = inv;
for c = 1:p.fieldcells
  inv{c} = 1 - 2*(a{c}.*conj(a{c}) > p.adapt);   %  1 for low amplitude
  a{c} = a{c}.^inv{c};                           % Invert field adaptively
end
a0 = a;                                          % initial cells
for iter = 1:p.iterations                        % Midpoint iteration loop
  for c = 1:p.fieldcells
     ac{c} = a{c}.^inv{c};
  end
  if p.setboundaries                            % if boundaries 
    [ac,p.boundval] = p.setbound(ac,p);         % set boundary values
  end                                           % end if boundaries
  for c = 1:p.fieldcells
     d = reshape(p.deriv{c}(ac{:},xi{:},p)*dt,p.d.ca{c});
     d = inv{c}.*d.*(a{c}.^(1-inv{c}));
     a{c} = a0{c} + d;
  end
end                                              % End iteration loop
for c = 1:p.fieldcells
  a{c} = 2*a{c}-a0{c};                           % Add derivatives
  a{c} = a{c}.^inv{c};                           % Invert field adaptively
end
p.t = p.t + dt;                                  % get endpoint time
a = p.prop(a,p);                                 % Propagate the fields
varargout = {1,2,2,1,1}; 
%  St. order = 1, det. order = 2, ipsteps = 2, vect = 1(Y), cell  = 1(Y)
end                                              % End function call    