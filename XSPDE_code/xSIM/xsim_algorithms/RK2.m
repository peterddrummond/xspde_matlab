function a  =  RK2(a,z,r)        
%   a = RK2(a,xi,dt,r)  propagates a step with second-order Runge-Kutta.   
%   Input: field a, noise xi, parameters r. Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

am  =  r.prop(a,r);                        %%propagates first RK2 field
dt  =  r.dtr;
a   =  r.prop(a + reshape(r.da(a,z,r)*dt,r.d.a),r);
r.t =  r.t + dt;                            %%Increment current time      
am  =  am + reshape(r.da(a,z,r)*dt,r.d.a);
a   =  0.5*(a + am); 
end                                         %%end function