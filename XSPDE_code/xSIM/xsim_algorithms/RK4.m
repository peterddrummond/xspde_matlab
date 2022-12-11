function a  =  RK4(a,xi,r)         
%   a = RK4(a,xi,r)  propagates a step with fourth-order Runge-Kutta.   
%   Input: field a, lattice r, noise xi.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

    dt = 0.5*r.dtr;
    ai = r.prop(a,r);
    d1 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%first derivative
    a1 = r.prop(a+d1,r);                           %%first estimate=ai+k1
    r.t = r.t+dt;                                  %%Increment time
    a2 = ai+reshape(r.da(a1,xi,r)*dt,r.d.a);       %%Second  estimate
    a  = a2;
    if r.setboundaries
        a  =  xsetbound(a,r);
    end
    d3 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%Third estimate
    r.t = r.t+dt;                                  %%Increment current time
    a =  r.prop(ai+2.*d3,r);                       %%Last field estimate
    d4 = reshape(r.da(a,xi,r)*dt,r.d.a);           %%fourth deriv
    a1 = r.prop((a1 + 2.*(a2 + d3))/3.,r);         %%Sum IP derivatives
    a = a1 + d4/3.;                                %%final algorithm
    if r.setboundaries
        a  =  xsetbound(a,r);
    end
end                                                %%end function
