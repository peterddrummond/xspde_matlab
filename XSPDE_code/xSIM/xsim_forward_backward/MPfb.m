function a  =  xMPfb(a,pr,xi,r)            
%   a = XMPFB(a,xi,r) propagates a forward-backward step using midpoints.   
%   Input: field a, previous iteration fields pr, noise xi, lattice r. 
%   Output: new field a. 
%   SPDE functions are licensed by Peter D. Drummond, (2020) - see License 

a = r.prop(a,r);                           %%Interaction picture propagator
p = 0.5*(pr{1}+pr{2});                     %%Previous iteration midpoint
dt = .5*r.dtr;                             %%Half time-step
r.t = r.t + dt;                            %%Increment current time
a0 = a;                                    %%Initialize iteration
for iter = 1:r.iterations                  %%Midpoint iteration loop
    d = reshape(r.da(a,p,xi,r)*dt,r.d.a);  %%Reshape midpoint derivative
    a = a0 + d;                            %%Increment fb field       
end                                        %%End iteration loop
a = a + d;                                 %%Increment fb field again
a = r.prop(a,r);                           %%Interaction picture propagator
end                                        %%end function 