function v = xgaussrandom(r)
%   v = XGAUSSRANDOM(r) generates a gaussian random matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates [randoms(1),nlattice] x-delta-correlated fields
%   Generates [randoms(2),nlattice] filtered k-delta-correlated fields
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
mrandoms = max(r.randoms(1),r.randoms(2));              %%Calculate max number
randoms = r.randoms(1)+r.randoms(2);
v = r.s.dx*randn(mrandoms,r.nlattice);                  %%Generate randoms
if r.randoms(2) > 0                                     %%Check if k-randoms
    kv = v(1:r.randoms(2),:);
    kv = reshape(kv,[r.randoms(2),r.d.int]);            %%reshape for fft
    for nd = 3:r.dimension+1                            %%loop over dimension
        kv = fft(kv,[],nd);                             %%inverse FFT
    end                                                 %%end if k-space noise  
    kv = kv.*r.infilt;                                  %%filter k-randoms 
    for nd = 3:r.dimension+1                            %%loop over dimension
        kv = ifft(kv,[],nd);                            %%inverse FFT
    end                                                 %%end dimension loop  
    kv = reshape(kv,r.randoms(2),r.nlattice);           %%reshape for XSPDE
    v(1+r.randoms(1):randoms,:) = kv;                   %%add to x-randoms
end                                                     %%End check k-randoms
end                                                     %%End function