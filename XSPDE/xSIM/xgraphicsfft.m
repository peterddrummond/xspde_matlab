function a  =  xgraphicsfft(a,trans,r)            
%   a = XGRAPHICSFFT(a,trans,r) selectively transforms spatial lattice fields.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'trans(i)' = 0 for space domain, = 1 for transform domain.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
a =reshape(a,r.d.fieldsplus);                   %%reshape to lattice
for nd = 2:r.dimension                          %%loop over space dimension
    if trans(nd-1) >0                           %%if FFT required
        a = fftshift(fft(a,[],2+nd)*r.kfact(nd),2+nd);%%Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a,r.d.aplus);                        %%reshape to flat array
end                                             %%end function