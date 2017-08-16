function a  =  xprop(a,r) 
%   a = XPROP(a,r) propagates a step in time using da. 
%   Input: field a, lattice r.
%   Output: new field a. 
%   Note, first two dimensions of a are components and ensembles.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
if r.propagator == 1                          %%no change required
    return                                    %%return field unaltered
else                                          %%change IS required
    a =reshape(a,r.d.fields);                 %%reshape to array
    dmax = r.dimension+2;                     %%maximum dimension
    for nd = 4:dmax                           %%loop over dimension
        a = fft(a,[],nd);                     %%take Fourier transform
    end                                       %%end loop over dimension
    a = r.propagator.*a;                      %%propagate in Fourier space
    for nd = 4:dmax                           %%loop over dimension
        a = ifft(a,[],nd);                    %%inverse Fourier transform
    end                                       %%end loop over dimension
    a = reshape(a, r.d.a);                    %reshape to matrix
end                                           %%end if
end                                           %%end xprop function