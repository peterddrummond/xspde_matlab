function a  =  xft(a,r) 
%   a = XFT(a,r) carries out a spatial inverse Fourier transform on a. 
%   Input: single component field a, lattice r.  Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
    a =reshape(a,r.d.int);                 %%reshape to array
    for nd = 2:r.dimension                 %%loop over dimension
        a = fft(a,[],nd);                  %%take Fourier transform
    end                                    %%end loop over dimensio
    a = reshape(a, [1,r.nlattice]);        %%reshape to matrix
end                                        %%end function