function a  =  xMPa(a,xi,r)            
%   a = XMPa(a,xi,r) propagates a step using the midpoint method 
%   includes anti-aliasing.   
%   Input: field a, noise xi, lattice r. 
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

    a(1:r.fields,:) = r.prop(a(1:r.fields,:),r);
    dt =  .5*r.dtr; 
    r.t = r.t + dt;                        %%Increment current time
    a0 = a;                                %%Initialize iteration
    for iter=1:r.iterations                %%Midpoint iteration loop
        d=reshape(r.da(a,xi,r)*dt,r.d.a);
        a(1:r.fields,:)=a0(1:r.fields,:)+d;
        if r.defines
            a = [a(1:r.fields,:);r.define(a,xi,r)]; %%Get defined fields
        end               
    end                                    %%End iteration loop
    
    a(1:r.fields,:) = r.prop(a(1:r.fields,:),r);%initial plus  derivative
    a(1:r.fields,:)=a(1:r.fields,:)+r.propa(d,r);
    a = [a(1:r.fields,:);r.define(a,xi,r)];
end                                        %%end function 