function w = xrandomgen(r)
%   xi = = XRANDOMGEN(r) generates random field matrix w
%   Input 'r' structure includes lattice dimensions
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
w = r.s.dx*randn(r.randoms(1),r.n.lattice);                 %%Calculate xnoise
if r.randoms(2) >0                                          %%check if k-noise 
     kw = r.s.dk*randn(r.randoms(2),r.n.lattice).*r.infilt; %% knoise
     kw = reshape(kw,[r.randoms(2),r.d.int]);               %%ireshape for fft 
     for nd = 3:r.dimension+1                               %%loop over dimension
         kw = ifft(kw,[],nd);                               %%inverse Fourier transform
     end                                                    %%end loop over dimension
     kw = reshape(kw,r.randoms(2),r.n.lattice);
     w = vertcat(w,kw);                                     %%add to x-noise
end                                                         %%End check if k-noise
end