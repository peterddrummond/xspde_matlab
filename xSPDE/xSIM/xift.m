function a  =  xift(a,r) 
%   a = XIFT(a,r) carries out a spatial inverse Fourier transform on a. 
%   Input: FT field a, input structure r.  Output: new field a, without FT.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
    sz = size(a);
    a =reshape(a,[sz(1),r.d.int]);         %%reshape to array 
    for nd = 2:r.dimension                 %%loop over dimension
        a = ifft(a,[],nd+1);               %%take Fourier transform
    end                                    %%end loop over dimensio
    a = reshape(a, [sz(1),r.nlattice]);    %%reshape to matrix
end                                        %%end  function