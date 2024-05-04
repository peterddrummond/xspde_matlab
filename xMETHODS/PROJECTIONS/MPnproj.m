function [a,varargout]  =  MPnproj(a,xi,p)            
%   a = MPnproj(a,xi,r) propagates a step using a midpoint method.
%   Includes a tangential derivative plus a final normal projection
%   Input: cell a, noise cell xi, parameters p
%   Outputs: cell a, code signatures varargout
%   Needs a geometry projection handle, p.project(d,a,n,(c,)p).
%   The projector function requires arguments of (d,a,n,(c,)p) where:
%   d is the derivative, c is the optional cell index
%   n=1 - returns a tangential projection of d near location a
%   n=2 - returns a normal projection of a{c} to the constraint  
%   Licensed by Peter D. Drummond, (2020) - see License


a =  p.prop(a,p);                                % Propagate the fields
dt  = p.dtr/2.0;                                 % Half time-step 
p.t = p.t + dt;                                  % get midpoint time
a0 = a; d = a;                                   % initial cells
for iter = 1:p.iterations                        % Midpoint iteration loop
  for c = 1:p.fieldcells
     d{c} = reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});
  end
  for c = 1:p.fieldcells
    a{c} = a0{c} + p.project(d{c},a,1,c,p);      % Tangential derivative;
  end                                            % End cell loop
end                                              % End iteration loop
for c = 1:p.fieldcells
  a0{c} = 2*a{c}-a0{c};                          % Add derivatives
end
for c = 1:p.fieldcells
  a{c} = p.project(0,a0,2,c,p);                  % Project normally
end
a = p.prop(a,p);                                 % Propagate the fields
varargout = {1,2,2,1,1};
%% Stoch. order = 1, det. order = 2, ipsteps = 2, vector = Y, cell  = Y
end                                              % End function call    