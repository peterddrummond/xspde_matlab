function a  =  xRK4(a,xi,r)         
%   a = XRK4(a,xi,r)  propagates a step with fourth-order Runge-Kutta.   
%   Input: field a, lattice r, noise xi.
%   Output: new field a. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

am = r.prop(a,r);                              %%Initialises midpoint RK4
d1 = reshape(r.da(a,xi,r)*r.dtr,r.d.a);        %%first derivative
d1 = 0.5*r.prop(d1,r);                         %%first derivative
r.t = r.t+0.5*r.dtr;                           %%Increment current time         
d2 = 0.5*reshape(r.da(am+d1,xi,r)*r.dtr,r.d.a);%%Second deriv
d3 = 0.5*reshape(r.da(am+d2,xi,r)*r.dtr,r.d.a);%%third deriv
r.t = r.t+0.5*r.dtr;                           %%Increment current time
a =  r.prop(am+2.*d3,r);                       %%Last field estimate
d4 = 0.5*reshape(r.da(a,xi,r)*r.dtr,r.d.a);    %%fourth deriv
d1 = (d1 + 2.*(d2 + d3))/3. ;                  %%Sum midpoint derivatives
a =  r.prop(am + d1,r) + d4/3.;                %%RK4IP final algorithm
end                                            %%end function