function a  =  xift(a,r) 
%   a = XIFT(a,r) carries out a spatial inverse Fourier transform on a. 
%   Input: FT field a, input structure r.  Output: new field a, without FT.
%   First dimension is the field index, last dimension is the ensemble
%   Note: definition is that a_out(x) = (1/N)sum(exp(ikx)a_in(k)) 
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License 
                                                            
for nd = 2:r.dimension                       %%loop over dimension
    a = ifft(a,[],nd);                       %%take Fourier transform
end                                          %%end loop over dimensio
end                                          %%end  function