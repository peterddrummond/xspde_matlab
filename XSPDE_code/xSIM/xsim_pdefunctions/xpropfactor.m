function kpropagator = xpropfactor (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Uses derivatives as either r.Dx or r.D{1}
%   Output: 'kpropagator' is a k-space multiplicative factor used by xprop.
%   For sde case uses a matrix propagator
%   First dimension is either the field index or 1
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

L =  r.linear(r);                                %%Get the linear term
if ~isequal(L, zeros(r.d.r))                     %%if L nonzero
  if r.dimension > 1                             %%If SPDE case
    sz = size(L);                                %%get size of L
    L = reshape(L,[sz(1),r.d.space]);            %%reshape to expected size
    kpropagator = exp(L*r.dt/(nc*r.ipsteps));    %%get propagation factor
  else                                           %%else sde
    kpropagator = expm(L*r.dt/(nc*r.ipsteps));   %%matrix propagator
  end                                            %%End if SPDE case
else                                             %%else L is zero
  kpropagator = 0;                               %%set 0 for no propagator
end                                              %%end if L nonzero
end
