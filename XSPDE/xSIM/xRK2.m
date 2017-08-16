function a  =  xRK2(a,z,r)        
%   a = XRK2(a,xi,dt,r)  propagates a step with second-order Runge-Kutta.   
%   Input: field a, noise xi, parameters r. Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

if r.defines
    am = r.prop(a(1:r.fields,:),r);         %%propagates first RK2 field 
    dt =  r.dtr;
    a(1:r.fields,:)=a(1:r.fields,:)+reshape(r.da(a,z,r)*dt,r.d.a);
    a(1:r.fields,:)=r.prop(a(1:r.fields,:),r);
    r.t=r.t+dt;                              %%Increment current time      
    a = [a(1:r.fields,:);r.define(a,z,r)];  %%Estimate field at last time
    am=am+reshape(r.da(a,z,r)*dt,r.d.a);
    a(1:r.fields,:) = 0.5*(a(1:r.fields,:)+am); 
    a = [a(1:r.fields,:);r.define(a,z,r)];
else
    am = r.prop(a,r);         %%propagates first RK2 field
    dt =  r.dtr;
    a=r.prop(a+reshape(r.da(a,z,r)*dt,r.d.a),r);
    r.t=r.t+dt;                              %%Increment current time      
    am=am+reshape(r.da(a,z,r)*dt,r.d.a);
    a = 0.5*(a+am); 
end
end                                        %%end function