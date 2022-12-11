function a  =  xEtn(a,xi,r)            
%   a = xEtn(a,xi,r) propagates a projected step using an Euler method.
%   Includes a tangential derivative and a final normal step
%   Input: field a, noise xi, data structure r 
%   Needs a geometry projection handle, projector r.project(d,a,n,r).
%   where:  n = 1 - tangential projection, n = 2 - normal projection 
%   Output: new field a, propagated and projected.
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License
d = r.da(a,xi,r)*r.dtr;                    %%Calculate estimated derivative
a = a + r.project(d,a,1,r) ;               %%Get tangential derivative
a = r.project(0,a,2,r) ;                   %%Project the field normally
end                                        %%End function call  