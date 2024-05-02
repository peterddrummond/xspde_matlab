function a  =  RK4J(a,xi,p)            
%   a = RK4J(a,xi,p) propagates a quantum jump step with RK4.   
%   Input: field a, noise w, parameters p. Output: new field a.
%   SDE functions are licensed by Peter D. Drummond, (2023) - see License 
dt  = 0.5*p.dtr;                                 %%Half time-step
a   = JumpB(a,p.dtr,p);                             %%Jump
da  = p.d.a;                                     %%Shape
a1  = a + dt*reshape(JumpA(a,xi,p),da);          %%first derivative
p.t = p.t+dt;                                    %%Increment time
a2  = a + dt*reshape(JumpA(a1,xi,p),da);         %%Second  estimate
a3  = a + 2*dt*reshape(JumpA(a2,xi,p),da);       %%Third estimate
p.t = p.t+dt;                                    %%Increment  time
a   = (a1 + 2.*a2 + a3 - a)/3.0;                 %%Sum three derivatives
a   = a+dt*reshape(JumpA(a3,xi,p),da)/3.0;       %%fourth deriv
%a   = JumpB(a,dt,p);                             %%Jump
if p.quantum == 1 
  a = a./sqrt(sum(conj(a).*a,1:p.nfields));      %%Project wavefunction
elseif p.quantum == 2
  a = a./trace(a);                               %%Project density
end                                              %%End projection
end                                              %%end function