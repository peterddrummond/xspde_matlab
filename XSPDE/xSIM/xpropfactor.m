function propagator = xpropfactor (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Uses derivatives as either r.Dx or r.D{1}
%   Output: 'propagator' is a lattice multiplicative factor used by xprop. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

L = r.linear(r);
sz = size(L);
if isequal(sz,r.d.a)
    if isequal(L,zeros(r.d.a))
        propagator = 1;
        return
    end
elseif sz(1) == r.d.a(1)   
    L = repmat(L(:,1),1,r.d.a(2));
else 
    L = repmat(L(1,1),r.d.a);
end
L = reshape(L,r.d.fields);  
propagator = exp(L*r.dt/(nc*r.ipsteps));  %%propagation factor
end