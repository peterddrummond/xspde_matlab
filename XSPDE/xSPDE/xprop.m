function a  =  xprop(a,r) 
%   a = XPROP(a,r) propagates a step in time using linear couplings. 
%   Input: field a, lattice r.
%   Output: new field a. 
%   Note, first two dimensions of a are components and ensembles.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 
                                                            
    a =reshape(a,r.d.ft);                         %%reshape to array
    dmax = r.dimension+1;                         %%maximum dimension
    for nd = 3:dmax                               %%loop over dimension
        a = fft(a,[],nd);                         %%take Fourier transform
    end                                           %%end loop over dimension
    a = r.propagator.*a;                          %%propagate in Fourier space
    for nd = 3:dmax                               %%loop over dimension
        a = ifft(a,[],nd);                        %%inverse Fourier transform
    end                                           %%end loop over dimension
    a = reshape(a, r.d.a);                        %%reshape to matrix
end                                               %%end xprop function