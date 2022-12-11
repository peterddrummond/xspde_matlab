function a  =  MP(a,xi,r)            
%   a = MP(a,xi,r) propagates a step using the midpoint method.   
%   Input: field a, noise xi, lattice r. 
%   Output: new field a.
%   Treats simultaneous equations with cell arrays of fields
%   SPDE functions are licensed by Peter D. Drummond, (2015) - see License 

    a = r.prop(a,r);                         %interaction picture propagate
    dt =  .5*r.dtr;                          %%half time-step 
    r.t = r.t + dt;                          %%Increment current time
    a0 = a;                                  %%Initialize iteration
    for iter = 1:r.iterations                %%Midpoint iteration loop
        d = reshape(r.da(a,xi,r)*dt,r.d.a);  %%Compute derivative
        a = a0+d;                            %%Compute midpoint
    end                                      %%End iteration loop
    a = a+d;                                 %%Compute endpoint
    a = r.prop(a,r);                         %interaction picture propagate
end                                          %%end function 