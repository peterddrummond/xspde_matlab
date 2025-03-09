function [ac,varargout]    =  RK2(ac,xi,p)        
%   ac = RK2(ac,xi,p) propagates a step with second-order Runge-Kutta + IP.   
%   Input:  cell 'ac',  noise 'xi', parameter structure 'p'.
%   Output: cell 'ac', including a derivative and a linear propagator. 
%   Licensed by Peter D. Drummond, (2024) - see License and manual.

a  = cell(1,p.fieldcells);
dt  = p.dtr;
for c = 1:p.fieldcells
    a{c} = ac{c} + reshape(p.deriv{c}(ac{:},xi{:},p)*dt,p.d.ca{c});      
end
p.t = p.t + dt;                                  % get endpoint time
ac =  p.prop(ac,p);                              % propagate to endpoint                                     
if p.setboundaries  
    [ac,p.boundval] = p.setbound(ac,p);
end  
a =  p.prop(a,p);                              % propagate to endpoint                                     
if p.setboundaries  
    [a,p.boundval] = p.setbound(a,p);
end  
for c = 1:p.fieldcells
    ac{c} = ac{c} + reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});
    ac{c} = (a{c} +  ac{c})/2;    
end 
varargout = {1,2,2,1,1};
%% Stoch. order = 1, det. order = 2, ipsteps = 1, vector = Y, cell  = Y
end                                              % end function