function a  =  Implicit(a,xi,r)            
%   a = Implicit(a,xi,r) propagates the implicit Ito-Euler method. 
%   Input: field a, noise xi, lattice r. 
%   Output: new field a. 
%   SPDE functions are licensed by Peter D. Drummond, (2015) - see License 

    a = r.prop(a,r);
    dt =  r.dtr; 
    r.t = r.t + dt;                        %%Increment current time
    a0 = a;                                %%Initialize iteration
    for iter=1:r.iterations                %%Implicit iteration loop 
        d=reshape(r.da(a,xi,r)*dt,r.d.a);
        a=a0+d;                    
    end                                    %%End iteration loop
end                                        %%end function 