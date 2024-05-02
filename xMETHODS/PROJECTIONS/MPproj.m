function [a,varargout]  =  MPproj(a,xi,p)            
%   a = MPproj(a,xi,r) propagates a step using a midpoint method.
%   Includes a tangential derivative without a final normal projection
%   Input: cell a, noise cell xi, data structure p
%   Outputs: cell a, code signatures varargout
%   Needs a geometry projection handle, p.project(d,a,n,r).
%   The projector function requires arguments of (d,a,n,r) where:
%   n=1 - returns a tangential projection of d near location a
%   Licensed by Peter D. Drummond, (2024) - see License


a =  p.prop(a,p);                                %%Propagate the fields
dt  = p.dtr/2.0;                                 %%Half time-step 
p.t = p.t + dt;                                  %%get midpoint time
a0 = a;                                          %%initial cells
for iter = 1:p.iterations                        %%Midpoint iteration loop
  for c = 1:p.fieldcells
     d = reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});
     d = p.project(d,a{c},1,p);                  %%Tangential derivative
     a{c} = a0{c} + d;
  end
end                                              %%End iteration loop
for c = 1:p.fieldcells
  a{c} = 2*a{c}-a0{c};                           %%Add derivatives
end
a = p.prop(a,p);                                 %%Propagate the fields
varargout = {1,2,2,1,1};
%% Stoch. order = 1, det. order = 2, ipsteps = 2, vector = Y, cell  = Y
end                                              %%End function call    