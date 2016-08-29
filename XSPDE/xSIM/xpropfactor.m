function propagator = xpropfactor (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Uses derivatives as either r.Dx or r.D{1}
%   Output: 'propagator' is a lattice multiplicative factor used by xprop. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

L = reshape(r.linear(r),r.d.fields);  
if isequal(L,zeros(r.d.a))
    propagator = 1;
else 
    propagator = exp(L*r.dt/(nc*r.ipsteps));  %%propagation factor
end
end
