function v = xgaussrandom(r)
%   v = XGAUSSRANDOM(r) generates a gaussian random matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates [randoms(1),nlattice] x-delta-correlated fields
%   plus [randoms(2),nlattice] (filtered) k-delta-correlated fields
%   Requires a filter function, r.rfilter(kv,r) in k-space.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
v = r.s.dx*randn([r.randoms(1),r.d.lattice]);       %%Generate randoms
if r.randoms(2) > 0                                 %%If k-randoms
    kv = randn([r.randoms(2),r.d.lattice]);         %%Normalize k-randoms
    kv = r.rfilter(kv,r);                           %%filter k-randoms
    kv = r.s.dk*xift(kv,r)/r.kfspace;               %%normalized IFFT
    v = [v;kv];                                     %%combined noises
end                                                 %%End if k-randoms
end                                                 %%End function