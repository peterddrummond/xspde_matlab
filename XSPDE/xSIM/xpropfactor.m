function propagator = xpropfactor (nc,r)                       
%   XPROPFACTOR(nc,r)  returns the interaction picture propagation factor.
%   Input:   check index 'nc', struct 'r' with coordinates.
%   Output: 'propagator' is a lattice multiplicative factor used by xprop. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 
  
if r.dimension < 2
    D=0;
else
    if r.numberaxis || r.dimension > 4
      D = cell(r.dimension-1);  
      for d=1:r.dimension-1
          D{d} = 1i*r.k{d};                            %%compute D, numeric
      end
    else
      D = struct('x',1i*r.kx,'y',1i*r.ky,'z',1i*r.kz); %%compute D, symbolic
    end
end
L = reshape(r.linear(D,r),r.d.fields);                    %%linear response
if isequal(L,zeros(r.d.fields))
    propagator = 1;
else 
    propagator = exp(L*r.dt/(nc*r.ipsteps));              %%propagation factor
end
end
