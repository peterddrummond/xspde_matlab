function [a,varargout]  =  Enproj(a,xi,p)            
%   a = Enproj(a,xi,r) propagates a projected step using an Euler method.
%   Includes a tangential derivative and a final normal step
%   Input: cell a, noise cell xi, data structure p
%   Outputs: cell a, code signatures varargout
%   Needs a geometry projection handle, p.project(d,a,n,r).
%   The projector function requires arguments of (d,a,n,r) where:
%   n=1 - returns a tangential projection of d near location a
%   n=2 - returns a normal projection of a to the constraint  
%   Licensed by Peter D. Drummond, (2020) - see License



dt  = p.dtr;                                     %%Time-step 
for c = 1:p.fieldcells
    d = reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});
    a{c} = a{c} + p.project(d,a{c},1,p);         %%Tangential derivative;
end                                              %%End iteration loop
p.t = p.t + dt;                                  %%get time
a = p.prop(a,p);                                 %%Propagate the fields
for c = 1:p.fieldcells
  a{c} = p.project(0, a{c} ,2,p);                %%Project normally
end
varargout = {1,1,1,1,1};
%% Stoch. order = 1, det. order = 1, ipsteps = 1, vector = Y, cell  = Y
end                                              %%End function call    