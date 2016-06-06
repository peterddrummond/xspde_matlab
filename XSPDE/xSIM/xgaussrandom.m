function v = xgaussrandom(r)
%   v = XGAUSSRANDOM(r) generates a gaussian random matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates [randoms(1),nlattice] x-delta-correlated fields
%   Generates [randoms(2),nlattice] filtered k-delta-correlated fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
v = r.s.dx*randn(r.randoms(1),r.nlattice);            %%Calculate x-randoms
if r.randoms(2) >0                                    %%Check if k-randoms 
     kv = r.s.dk*randn(r.randoms(2),r.nlattice);      %%generate k-randoms
     kv = kv.*r.infilt;                               %%filter k-randoms 
     kv = reshape(kv,[r.randoms(2),r.d.int]);         %%reshape for FFT 
     for nd = 4:r.dimension+2                         %%loop over dimension
         kv = ifft(kv,[],nd);                         %%inverse FFT
     end                                              %%end dimension loop  
     kv = reshape(kv,r.randoms(2),r.nlattice);        %%reshape for XSPDE
     v = vertcat(v,kv);                               %%add to x-randoms
end                                                   %%End check k-randoms
end                                                   %%End function