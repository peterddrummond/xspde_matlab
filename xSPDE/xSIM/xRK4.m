function ao  =  xRK4(a,xi,r)         
%   a = XRK4(a,xi,r)  propagates a step with fourth-order Runge-Kutta.   
%   Input: field a, lattice r, noise xi.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

    dt = 0.5*r.dtr;
    a1 = a(1:r.fields,:); 
    am = r.prop(a1,r);
    d1 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%first derivative
    a1 = r.prop(a1+d1,r);                          %%first estimate
    r.t = r.t+dt;                                  %%Increment time
    a = [a1;r.define(a,xi,r)];
    a2 = am+reshape(r.da(a,xi,r)*dt,r.d.a);        %%Second  estimate
    a(1:r.fields,:)=a2;
    if r.setboundaries
        a  =  xsetbound(a,r);
    end
    a = [a(1:r.fields,:);r.define(a,xi,r)];
    d3 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%Third estimate
    r.t = r.t+dt;                                  %%Increment current time
    a(1:r.fields,:) =  r.prop(am+2.*d3,r);         %%Last field estimate
    a = [a(1:r.fields,:);r.define(a,xi,r)];
    d4 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%fourth deriv
    a1 = r.prop((a1 + 2.*(a2 + d3))/3.,r);         %%Sum IP derivatives
    a(1:r.fields,:) = a1 + d4/3.;                  %%final algorithm
    if r.setboundaries
        a  =  xsetbound(a,r);
    end
    ao = [a(1:r.fields,:);r.define(a,xi,r)];
end                                                %%end function
