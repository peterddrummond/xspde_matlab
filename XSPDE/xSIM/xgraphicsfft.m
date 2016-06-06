function a  =  xgraphicsfft(a,r,tr)            
%   a = XGRAPHICSFFT(a,r,tr) selectively transforms spatial lattice fields.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'tr(i)' = 0 for space domain, = 1 for transform domain.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
a =reshape(a, r.d.fields);                          %%reshape to lattice
dmax = ndims(a);                                %%get a dimension
for nd = 3:dmax                                 %%loop over dimension
    if tr(nd-2) >0                              %%if FFT required
        a = fft(a,[],nd)*r.kfact(nd-1);         %%take Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a, r.d.a);                           %%reshape to flat array
end                                             %%end function