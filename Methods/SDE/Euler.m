function [ac,varargout]  =  Euler(ac,xi,p)         
%   ac = Euler(a,xi,p) propagates a step using the Euler + IP method.   
%   Input:  cell of fields 'ac',  noises 'xi', parameter structure 'p'.
%   Output: cell of fields 'ac',  with derivative and linear propagator.
%   Licensed by Peter D. Drummond, (2024) - see License and manual.
for c = 1:p.fieldcells
    ac{c} = ac{c} + reshape(p.deriv{c}(ac{:},xi{:},p)*p.dtr,p.d.ca{c});      
end
ac = p.prop(ac,p);                               %%Propagate the fields
varargout = {1,1,1,1,1}; 
%% St. order = 1, det. order = 1, ipsteps = 1, vect = 1(Y), cell  = 1(Y)
end                                              %%End function call                                                      