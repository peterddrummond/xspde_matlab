function kpropagatora = xpropfactora (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   includes anti-aliasing
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
            kpropagatora = 0;              %%Set default propagator = 0 
        else
            kpropagatora = 1;              %%Set default propagator = 1
        end
        return
    end
elseif sz(1) ~= r.d.a(1) || sz(2) ~= 1
    xfcheck ('linear',0,L,r.d.a);         %%Error message for wrong size
else  
    L = repmat(L(:,1),1,r.d.a(2));        %%Repeat over remaining indices
end
L = reshape(L,r.d.fields);
al=1;
for d=2:r.dimension
    al = al*sinc(0.5*r.k{d}*r.dx(d));
end
if r.broadcast == 1                       %%Check the broadcast option
    L1 = L(:,1,:);                        %%Use same factor in ensemble   
    L = reshape(L1,r.d.kfields);
end
kpropagatora = al*exp(L*r.dt/(nc*r.ipsteps));  %%propagation factor
end