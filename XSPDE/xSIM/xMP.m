function a  =  xMP(a,xi,r)            
%   a = XMP(a,xi,r) propagates a step using the midpoint method.   
%   Input: field a, lattice r, noise xi.
%   Output: new field a. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

if  r.defines
    a0 = a(1:r.fields,:);
    if ~isequal(r.propagator,1)
        a0 = r.prop(a0,r);
    end
    r.t = r.t + .5*r.dtr;                        %%Increment current time
    amid = a0;                                   %%Initialize iteration
    for iter=1:r.iterations                      %%Midpoint iteration loop
        a = [amid;r.define(a,xi,r)];
        d1 = r.da(a,xi,r)*r.dtr;                 %%get derivative fields
        d1 = reshape(d1,r.d.a);                  %%reshape fields  
        amid = a0 + 0.5*d1;                      %%midpoint fields
    end
    amid = a0 + d1;                              %%End iteration loop
    if ~isequal(r.propagator,1)
        amid = r.prop(amid,r);
    end
    r.t = r.t + .5*r.dtr;                        %%Increment current time
    a = [amid;r.define(a,xi,r)];
else
    if ~isequal(r.propagator,1)
        a = r.prop(a,r);
    end
    r.t = r.t + .5*r.dtr;                        %%Increment current time
    a0 = a;                                      %%Initialize iteration
    for iter=1:r.iterations                      %%Midpoint iteration loop   
        d1 = r.da(a,xi,r)*r.dtr;                 %%get derivative fields
        d1 = reshape(d1,r.d.a);                  %%reshape fields  
        a = a0 + 0.5*d1;                         %%midpoint fields                     
    end                                          %%End iteration loop
    a = a0 + d1;                                 %%initial plus  derivative
    if ~isequal(r.propagator,1)
        a = r.prop(a,r);
    end
end                                              %%end if   
end                                              %%end function
