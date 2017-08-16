function ao  =  xRK4(a,xi,r)         
%   a = XRK4(a,xi,r)  propagates a step with fourth-order Runge-Kutta.   
%   Input: field a, lattice r, noise xi.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

dt = 0.5*r.dtr;
if  ~r.defines && isequal(r.propagator,1)          %%simplest case
    d1 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%first derivative
    r.t = r.t+dt;                                  %%Increment current time         
    d2 = reshape(r.da(a+d1,xi,r)*dt,r.d.a);        %%Second deriv
    d3 = reshape(r.da(a+d2,xi,r)*dt,r.d.a);        %%third deriv
    r.t = r.t+dt;                                  %%Increment current time
    ao =  a+2.*d3;                                 %%Last field estimate
    d4 = reshape(r.da(ao,xi,r)*dt,r.d.a);          %%fourth deriv
    ao =   a + (d1 + 2.*(d2 + d3)+d4)/3. ;         %%Sum  derivatives
    fprintf('a=%12.8f\n',a(1,1));
    fprintf('d1=%12.8f\n',2*d1(1,1));
    fprintf('d2=%12.8f\n',2*d2(1,1));
    fprintf('d3=%12.8f\n',2*d3(1,1));
    fprintf('d4=%12.8f\n',2*d4(1,1));
    fprintf('aout=%12.8f\n',ao(1,1));
else
    am = r.prop(a(1:r.fields,:),r);
    d1 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%first derivative
    d1 = r.prop(d1,r);                             %%first derivative
    r.t = r.t+dt;                                  %%Increment current time
    a(1:r.fields,:)=am+d1;
    a = [a(1:r.fields,:);r.define(a,xi,r)];
    d2 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%Second deriv
    a(1:r.fields,:)=am+d2;
    a = [a(1:r.fields,:);r.define(a,xi,r)];
    d3 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%third deriv
    r.t = r.t+dt;                                  %%Increment current time
    a(1:r.fields,:) =  r.prop(am+2.*d3,r);         %%Last field estimate
    a = [a(1:r.fields,:);r.define(a,xi,r)];
    d4 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%fourth deriv
    d1 = (d1 + 2.*(d2 + d3))/3. ;                  %%Sum IP derivatives
    a(1:r.fields,:) = r.prop(am + d1,r) + d4/3.;   %final algorithm
    ao = [a(1:r.fields,:);r.define(a,xi,r)];
end
end                                                %%end function
