
function z = xphasenoise(r)
%   z = = XPHASENOISE(r) generates random-phase noise matrix z
%   Input 'r' structure includes lattice dimensions
%   Generates [noises(1),nlattice] x-t-delta-correlated fields
%   Generates [noises(2),nlattice] k-t-delta-correlated filtered fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
z = r.s.dxt*exp(2*pi*1i*rand(r.noises(1),r.nlattice)); %%Calculate  xnoise
if r.noises(2) > 0                                     %%if k-space noise
    zk = r.s.dkt*exp(2*pi*1i*rand(r.noises(1),r.nlattice));%%knoise generate
    zk = zk.*r.noisefilt;                               %%knoise filter
    zk = reshape(zk,[r.noises(2),r.d.int]);             %%reshape for fft 
    for nd = 3:r.dimension+1                            %%loop over dimension
        zk = ifft(zk,[],nd);                            %%inverse FFT
    end                                                 %%end if k-space noise  
    zk = reshape(zk,r.noises(2),r.nlattice);            %%reshape back
    z = vertcat(z,zk);                                  %%add to x-noise
end                                                     %%end if k-space noise  
end