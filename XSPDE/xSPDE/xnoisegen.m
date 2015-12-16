
function z = xnoisegen(r)
%   z = = XNOISEGEN(r) generates random noise matrix z
%   Input 'r' structure includes lattice dimensions
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
z = r.s.dxt*randn(r.noises(1),r.n.lattice);             %%Calculate  xnoise
if r.noises(2) >0                                       %%if k-space noise
    zk = r.s.dkt*randn(r.noises(2),r.n.lattice);        %%knoise
    zk = zk.*r.noisefilt;                               %%knoise filter
    zk = reshape(zk,[r.noises(2),r.d.int]);             %%reshape for fft 
    for nd = 3:r.dimension+1                            %%loop over dimension
        zk = ifft(zk,[],nd);                            %%inverse Fourier transform
    end                                                 %%end if k-space noise  
    zk = reshape(zk,r.noises(2),r.n.lattice);           %%reshape back
    z = vertcat(z,zk);                                  %%add to x-noise
end                                                     %%end if k-space noise  
end