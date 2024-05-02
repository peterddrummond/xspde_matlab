function [a,varargout]   =  Implicit(a,xi,p)            
%   a = Implicit(ac,xi,p) propagates a step using the Implicit + IP method.   
%   Input:  cell of fields 'ac',  noises 'xi', parameter structure 'p'.
%   Output: cell of fields 'ac',  with derivative and linear propagator. 
%   Licensed by Peter D. Drummond, (2024) - see License and manual.

a =  p.prop(a,p);                                %%Propagate the fields
dt  = p.dtr;                                     %%time-step 
p.t = p.t + dt;                                  %%get final time
a0 = a;                                          %%initial cells
for iter = 1:p.iterations                        %%Midpoint iteration loop
  for c = 1:p.fieldcells
     a{c} = a0{c} + reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});      
  end
end                                              %%End iteration loop
varargout = {1,1,1,1,1}; 
%% St. order = 1, det. order = 1, ipsteps = 1, vect = 1(Y), cell  = 1(Y)
end                                              %%End function call    