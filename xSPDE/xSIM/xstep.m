function a  =  xstep(a,aderiv,dt,xi,r) 
%   a = XSTEP(a,aderiv,r) propagates a step in time using da. 
%   Input: starting point a.
%   Evaluate deriv at aderiv, step dt, noise xi, lattice r.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
d = reshape(r.da(aderiv,xi,r)*dt,r.d.a);     %%Get deriv
a(1:r.fields,:) = a(1:r.fields,:)+d;
if r.defines
      a = [a(1:r.fields,:);r.define(a,xi,r)];%%Get defined fields
end
end                                           %%end xprop function