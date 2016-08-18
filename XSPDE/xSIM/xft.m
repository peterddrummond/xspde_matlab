function a  =  xft(a,r) 
%   a = XFT(a,r) carries out a spatial Fourier transform on a. 
%   INput: field a, input structure r.  Output: new field a, with FFT. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
    sz = size(a);
    a =reshape(a,[sz(1),r.d.int]);         %%reshape to array
    for nd = 2:r.dimension                 %%loop over dimension
        a = fft(a,[],nd+2);                %%take Fourier transform
    end                                    %%end loop over dimensio
    a = reshape(a, [sz(1),r.nlattice]);    %%reshape to matrix
end                                        %%end function