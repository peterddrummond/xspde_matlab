function z = xgaussnoise(r)
%   z = XGAUSSNOISE(r) generates Gaussian noise matrix z
%   Input 'r' structure includes lattice dimensions
%   Generates [noises(1),nlattice] x-delta-correlated noise fields
%   plus [noises(2),nlattice]  k-delta-correlated noise fields
%   Requires a filter function, r.nfilter(kv,r) in k-space.
%   Correlated noise fields generated by k-filtering plus inverse FFT
%   Output noises transformed to x-space, and used in stochastic equations
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

z = r.s.dxt*randn([r.noises(1),r.d.lattice]);           %%Get random gauss
if r.noises(2) > 0                                      %%if k-space noise
    zk = randn([r.noises(2),r.d.lattice]);              %%scale knoise
    zk = r.nfilter(zk,r);                               %%knoise filter
    zk = r.s.dkt*xift(zk,r)/r.kfspace;                  %%inverse FFT
    z=[z;zk];                                           %%add to noise
end                                                     %%end if k-noise
end