function a  =  xgft(a,trans,r)            
%   a = XGFT(a,trans,r) selectively transforms spatial lattice fields,
%   using conventional normalization and momentum component ordering.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'trans(i)' = 0 for no transform = 1 for forward transform.
%   First dimension is the field index, last dimension is the ensemble
%   OBSOLETE - not currently used in xspde
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
a =reshape(a,r.d.fieldsplus);                   %%reshape to lattice
for nd = 2:r.dimension                          %%loop over space dimension
    if trans(nd) > 0                            %%if FFT required
        a = fftshift(fft(a,[],1+nd)*r.kfact(nd),1+nd);%%Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a,r.d.aplus);                        %%reshape to flat array
end                                             %%end function