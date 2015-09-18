function a  =  xMP(a,xi,dt,r)            
%   a = XMP(a,xi,dt,r) propagates a step using the midpoint method.   
%   Input: field a, lattice r, noise xi, step dt.
%   Output: new field a. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

a = r.prop(a,r);
r.t = r.t + .5*dt;                            %%Increment current time
amid = a;                                     %%Initialize iteration
for iter=1:r.iterations                       %%Midpoint iteration loop   
     d1 = r.da(a,xi,r)*dt;                    %%get derivative fields
     d1 = reshape(d1,r.d.a);                  %%reshape fields  
     a = amid + 0.5*d1;                       %%midpoint fields                     
end                                           %%End iteration loop
a = amid + d1;                                %%initial plus  derivative
a = r.prop(a,r);                              %%Propagate midpoint field
end                                           %%end function
