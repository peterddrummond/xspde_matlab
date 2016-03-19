function v = xrandomphase(r)
%   v = XRANDOMPHASE(r) generates random phase matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates [randoms(1),nlattice] x-delta-correlated fields
%   Generates [randoms(2),nlattice] filtered k-delta-correlated fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
v = r.s.dx*exp(2*pi*1i*rand(r.noises(1),r.nlattice)) %%Calculate x-randoms
if r.randoms(2) >0                                    %%Check if k-randoms 
     kv = r.s.dk*exp(2*pi*1i*rand(r.noises(1),r.nlattice));%%gen k-randoms
     kv = kv.*r.infilt;                               %%filter k-randoms 
     kv = reshape(kv,[r.randoms(2),r.d.int]);         %%reshape for FFT 
     for nd = 3:r.dimension+1                         %%loop over dimension
         kv = ifft(kv,[],nd);                         %%inverse FFT
     end                                              %%end dimension loop  
     kv = reshape(kv,r.randoms(2),r.nlattice);        %%reshape for XSPDE
     v = vertcat(v,kv);                               %%add to x-randoms
end                                                   %%End check k-randoms
end                                                   %%End function