function kpropagator = xpropfactor (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Uses derivatives as either r.Dx or r.D{1}
%   Output: 'kpropagator' is a k-space multiplicative factor used by xprop.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

L = r.linear(r);
sz = size(L);                             %%Get size of linear propagator
if isequal(sz,r.d.a)                      %%Check for input matching fields
    if isequal(L,zeros(r.d.a))            %%Check for default L = zero
        if ~r.setboundaries 
            kpropagator = 0;              %%Set default propagator = 0 
        else
            kpropagator = 1;              %%Set default propagator = 1
        end
        return
    end
elseif sz(1) ~= r.d.a(1) || sz(2) ~= 1
    xfcheck ('linear',0,L,r.d.a);         %%Error message for wrong size
else  
    L = repmat(L(:,1),1,r.d.a(2));        %%Repeat over remaining indices
end
L = reshape(L,r.d.fields); 
if r.broadcast == 1                       %%Check the broadcast option
    L1 = L(:,1,:);                        %%Use same factor in ensemble   
    L = reshape(L1,r.d.kfields);
end
kpropagator = exp(L*r.dt/(nc*r.ipsteps));  %%propagation factor
end