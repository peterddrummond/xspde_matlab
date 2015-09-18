function propagator = xpropfactor (nc,r)                       
%   r = XPROPFACTOR(nc,r)  sets the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Output: 'propagator' is a lattice multiplicative factor used by xprop. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

  D.x=0;                                            %%initialize D 
  if r.dimension > 1                                %%if space dimensions 
    D = struct('x',1i*r.kx,'y',1i*r.ky,'z',1i*r.kz);%%compute D
  end                                               %%end if space 
  L = reshape(r.linear(D,r),r.d.ft);                %%linear response
  propagator = exp(L*r.dt/(nc*r.ipsteps));          %%propagation factor
end
