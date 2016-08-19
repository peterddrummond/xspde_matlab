function a  =  xRK2(a,xi,r)        
%   a = XRK2(a,xi,dt,r)  propagates a step with second-order Runge-Kutta.   
%   Input: field a, noise xi, parameters r.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

am = r.prop(a,r);                           %%in.propagates first RK2 field

a = a+reshape(r.da(a,xi,r)*r.dtr,r.d.a);    %%First RK2 field
a = r.prop(a,r);                            %%in.propagates field + deriv
r.t=r.t+r.dtr;                              %%Increment current time      
am = am+reshape(r.da(a,xi,r)*r.dtr,r.d.a);  %%Second RK2 
a = 0.5*(a+am);                             %%average of derivatives
end                                         %%end function