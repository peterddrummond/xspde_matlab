function v = xrandomgen(r)
%   v = XRANDOMGEN(r) generates random field matrix v
%   Input 'r' structure includes lattice dimensions
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
v = r.s.dx*randn(r.randoms(1),r.n.lattice);                 %%Calculate xnoise
if r.randoms(2) >0                                          %%check if k-noise 
     kv = r.s.dk*randn(r.randoms(2),r.n.lattice).*r.infilt; %% knoise
     kv = reshape(kv,[r.randoms(2),r.d.int]);               %%ireshape for fft 
     for nd = 3:r.dimension+1                               %%loop over dimension
         kv = ifft(kv,[],nd);                               %%inverse Fourier transform
     end                                                    %%end loop over dimension
     kv = reshape(kv,r.randoms(2),r.n.lattice);
     v = vertcat(v,kv);                                     %%add to x-noise
end                                                         %%End check if k-noise
end