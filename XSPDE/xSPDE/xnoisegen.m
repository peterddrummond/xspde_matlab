function xi = xnoisegen(r)
%   xi = = XNOISEGEN(r) generates random noise matrix xi
%   Input 'r' structure includes lattice dimensions
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
xi = r.s.dxt*randn(r.noises(1),r.n.lattice);            %%Calculate  xnoise
if r.noises(2) >0                                       %%if k-space noise
    xik = r.s.dkt*randn(r.noises(2),r.n.lattice);       %%knoise
    xik = xik.*r.noisefilt;                             %%knoise filter
    xik = reshape(xik,[r.noises(2),r.d.int]);           %%reshape for fft 
    for nd = 3:r.dimension+1                            %%loop over dimension
        xik = ifft(xik,[],nd);                          %%inverse Fourier transform
    end                                                 %%end if k-space noise  
    xik = reshape(xik,r.noises(2),r.n.lattice);         %%reshape back
    xi = vertcat(xi,xik);                               %%add to x-noise
end                                                     %%end if k-space noise  
end