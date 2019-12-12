function v = xgaussrandom(r)
%   v = XGAUSSRANDOM(r) generates a gaussian random matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates [randoms(1),nlattice] x-delta-correlated fields
%   plus [randoms(2),nlattice] (filtered) k-delta-correlated fields
%   Requires a filter function, r.rfilter(kv,r) in k-space.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
randoms = r.randoms(1)+r.randoms(2);                %%Calculate max number
v = randn(randoms,r.nlattice);                      %%Generate randoms
v(1:r.randoms(1),:)=r.s.dx*v(1:r.randoms(1),:);     %%Normalize x-randoms
if r.randoms(2) > 0                                 %%If k-randoms
    kv = r.s.dk*v(r.randoms(1)+1:randoms,:);        %%Normalize k-randoms
    kv = r.rfilter(kv,r);                           %%filter k-randoms
    kv = xift(kv,r)/r.kfspace;                      %%normalized IFFT
    v(1+r.randoms(1):randoms,:) = kv;               %%add to x-randoms    
end                                                 %%End if k-randoms
end                                                 %%End function