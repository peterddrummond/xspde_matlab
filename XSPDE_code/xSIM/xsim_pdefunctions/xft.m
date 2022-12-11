function a  =  xft(a,r) 
%   a = XFT(a,r) carries out a spatial Fourier transform on a. 
%   Input: field a, input structure r.  Output: new field a, with FFT.
%   First dimension is the field index, last dimension is the ensemble
%   Note: definition is that a_out(x) = sum(exp(-ikx)a_in(k)) 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
a = xshape(a,r);                           %%reshape to array
for nd = 2:r.dimension                     %%loop over dimension
    a = fft(a,[],nd);                      %%take Fourier transform
end                                        %%end loop over dimensio
end                                        %%end function