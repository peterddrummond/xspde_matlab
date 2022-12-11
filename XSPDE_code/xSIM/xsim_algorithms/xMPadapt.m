function a  =  xMPadapt(a,xi,r)            
%   a = XMPadapt(a,xi,r) propagates adaptive steps using a midpoint method.   
%   Input: field a, noise xi, lattice r. 
%   Output: new field a.
%   Integrates large amplitudes adaptively using an inverse method
%   Requires "r.adapt" to indicate critical |a|^2 to switch to inverse
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = r.prop(a,r);                           %%linear propagation
inv = 1 - 2*(a.*conj(a) > r.adapt);        %% 1 for low amplitude, -1 large
a = a.^inv;                                %%Invert field adaptively
dt = .5*r.dtr;                             %%Half time-step
r.t = r.t + dt;                            %%Increment current time
a0 = a;                                    %%Initialize iteration
for iter = 1:r.iterations                  %%Midpoint iteration loop
 	d = r.da(a.^inv,xi,r)*dt;              %%Usual derivative, non-adapted
    d = inv.*d.*(a.^(1-inv));              %%Adapted derivative
    a=a0+d;                                %%Next midpoint in iteration
end                                        %%End iteration loop
a=a+d;                                     %%Full step
a = a.^inv;                                %%Re-invert field adaptively
a = r.prop(a,r);                           %%Linear propagation
end                                        %%End function 