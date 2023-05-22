function a  =  MPnproj(a,xi,r)            
%   a = MPnproj(a,xi,r) propagates a step using a midpoint method.
%   Includes a tangential derivative plus a final normal projection
%   Input: field a, noise xi, data structure r 
%   Needs a geometry projection handle, r.project(d,a,n,r).
%   The projector function requires arguments of (d,a,n,r) where:
%   n=1 - returns a tangential projection of d near location a
%   n=2 - returns a normal projection of a to the constraint  
%   Output: new field a, propagated and projected.
%   SPDE functions are licensed by Peter D. Drummond, (2020) - see License

dt =  r.dtr/2;                             %%Half time-step for iteration
r.t = r.t + dt;                            %%Increment to get midpoint time
a0 = a;                                    %%Initialize iteration
for iter = 1:r.iterations                  %%Midpoint iteration loop
    d = r.da(a,xi,r)*dt;                   %%Calculate estimated derivative
    d = r.project(d,a,1,r);                %%Get tangential derivative
    a = a0 + d;                            %%Get the next midpoint estimate
end                                        %%End iteration loop
a = a + d;                                 %%Get the endpoint field estimate
a = r.project(0,a,2,r) ;                   %%Project the field normally
end                                        %%End function call  