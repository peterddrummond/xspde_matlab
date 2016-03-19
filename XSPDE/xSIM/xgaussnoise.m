
function z = xnoisegen(r)
%   z = = XNOISEGEN(r) generates random noise matrix z
%   Input 'r' structure includes lattice dimensions
%   Generates [noises(1),nlattice] delta-correlated noise fields
%   Generates [noises(2),nlattice] filtered noise  fields
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
z = r.s.dxt*randn(r.noises(1),r.nlattice);              %%Calculate  xnoise
if r.noises(2) > 0                                      %%if k-space noise
    zk = r.s.dkt*randn(r.noises(2),r.nlattice);         %%knoise generate
    zk = zk.*r.noisefilt;                               %%knoise filter
    zk = reshape(zk,[r.noises(2),r.d.int]);             %%reshape for fft 
    for nd = 3:r.dimension+1                            %%loop over dimension
        zk = ifft(zk,[],nd);                            %%inverse FFT
    end                                                 %%end if k-space noise  
    zk = reshape(zk,r.noises(2),r.nlattice);            %%reshape back
    z = vertcat(z,zk);                                  %%add to x-noise
end                                                     %%end if k-space noise  
end