function sz  =  xsize(o,r)                 
%   a  =  XSIZE(o,r) reports xSPDE array sizes with trailing ones included
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sz  = size(o);                               %%initial size
sz(end+1:r.dimension+1+r.indext)  = 1;       %%add trailing ones
end                                          %%end function