function v = xrandomgen(r)
%   v = XRANDOMGEN(r) generates random field matrix v
%   Input 'r' structure includes lattice dimensions
%   Generates r.randoms(1) delta-correlated random fields
%   Generates r.randoms(2) filtered random fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
v = r.s.dx*randn(r.randoms(1),r.nlattice);            %%Calculate x-randoms
if r.randoms(2) >0                                    %%Check if k-randoms 
     kv = r.s.dk*randn(r.randoms(2),r.nlattice);      %%generate k-randoms
     kv = kv.*r.infilt;                               %%filter k-randoms 
     kv = reshape(kv,[r.randoms(2),r.d.int]);         %%reshape for FFT 
     for nd = 3:r.dimension+1                         %%loop over dimension
         kv = ifft(kv,[],nd);                         %%inverse FFT
     end                                              %%end dimension loop  
     kv = reshape(kv,r.randoms(2),r.nlattice);        %%reshape for XSPDE
     v = vertcat(v,kv);                               %%add to x-randoms
end                                                   %%End check k-randoms
end                                                   %%End function