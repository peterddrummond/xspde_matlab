function a  =  xRK2(a,xi,dt,r)        
%   a = XRK2(a,xi,dt,r)  propagates a step with second-order Runge-Kutta.   
%   Input: field a, lattice r, noise xi, step dt.
%   Output: new field a. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

am = in.prop(a,r);                          %%in.propagates first RK2 field
a = a+reshape(in.da(a,xi,r)*dt,r.d.a);      %%First RK2 field
a = in.prop(a,ft);                          %%in.propagates field + deriv
r.t=r.t+dt;                                 %%Increment current time      
am = am+reshape(in.da(a,xi,r)*dt,r.d.a);    %%Second RK2 
a = 0.5*(a+am);                             %%average of derivatives
end                                         %%end function
