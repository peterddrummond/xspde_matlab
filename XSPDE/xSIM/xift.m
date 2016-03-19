function a  =  xift(a,r) 
%   a = XIFT(a,r) carries out a spatial inverse Fourier transform on a. 
%   Input: single component field field a, lattice r.  Output: inverse ft a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
    a =reshape(a,r.d.int);                 %%reshape to array
    for nd = 2:r.dimension                 %%loop over dimension
        a = ifft(a,[],nd);                 %%take Fourier transform
    end                                    %%end loop over dimensio
    a = reshape(a, [1,r.nlattice]);        %%reshape to matrix
end                                        %%end  function