function [a,varargout]  =  RK4(a,xi,p)         
%   a = RK4(a,xi,p)  propagates a step with fourth-order Runge-Kutta IP.   
%   Input: cell array of fields a, cell array of noises xi, parameters p.
%   Output: new cell array of fields a. 
%   Licensed by Peter D. Drummond, (2024) - see License

dt  = 0.5*p.dtr;
a1  = cell(1,p.fieldcells); a2 = a1; a3 = a1;  %%Initialize cells
for c = 1:p.fieldcells
  a1{c} = a{c}+reshape(p.deriv{c}(a{:},xi{:},p)*dt,p.d.ca{c});
end
a1  = p.prop(a1,p);                              %%First estimate
a   = p.prop(a,p);
p.t = p.t+dt;                                    %%Increment time
for c = 1:p.fieldcells
  a2{c} = a{c}+reshape(p.deriv{c}(a1{:},xi{:},p)*dt,p.d.ca{c});
end
if p.setboundaries
  a2 = p.setbound(a2,p);
end
for c = 1:p.fieldcells
  a3{c} = a{c}+reshape(p.deriv{c}(a2{:},xi{:},p)*2*dt,p.d.ca{c});
  a{c} = (a1{c} + 2*a2{c} + a3{c} - a{c})/3;
end
p.t = p.t+dt;                                    %%Increment time
a   = p.prop(a,p);
a3  = p.prop(a3,p);                              %%Last field estimate
for c = 1:p.fieldcells
  a{c} = a{c}+reshape(p.deriv{c}(a3{:},xi{:},p)*dt/3,p.d.ca{c});
end
if p.setboundaries
  a = p.setbound(a,p);
end
varargout = {1,4,2,1,1};
%% Stoch. order = 1, det. order = 4, ipsteps = 2, vector = Y, cell  = Y                             
end                                              %%end function