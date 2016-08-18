function a  =  xift(a,r) 
%   a = XIFT(a,r) carries out a spatial inverse Fourier transform on a. 
%   INput: FT field a, input structure r.  Output: new field a, without FT. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
    sz = size(a);
    a =reshape(a,[sz(1),r.d.int]);         %%reshape to array 
    for nd = 2:r.dimension                 %%loop over dimension
        a = ifft(a,[],nd+2);               %%take Fourier transform
    end                                    %%end loop over dimensio
    a = reshape(a, [sz(1),r.nlattice]);    %%reshape to matrix
end                                        %%end  function