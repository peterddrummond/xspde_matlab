function o  =  xshape(o,indext,r)                 
%   a  =  XSHAPE(o,r) shapes a field to standard xSPDE dimensions. 
%   Input: field variable 'o', lattice 'r'.
%   Checks if indext is set for possible time-index
%   Can shape either xSPDE fields or averages
%   Broadcasts over an ensemble if needed
%   Leaves fields with non-standard total components unchanged.
%   First dimension is the field or line index 
%   if r.indext is zero, the last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sz  = size(o);                               %%initial size
sz2 = sz(2:end);                             %%size of indices > 1
if indext                                    %%if time index is set
    lattice = r.points;                      %%lattice is space-time
else                                         %%if time index is not set
    lattice = r.d.lattice;                   %%lattice is space-ensemble
end                                          %%end if time index is set
if isequal(sz2,lattice)                      %%if lattice OK
    return                                   %%return
end                                          %%end if lattice OK
sz2 = prod(sz2);                             %%store remaining components
lsize = prod(lattice);                       %%total lattice size
if ~isequal(sz2,lsize)                       %%if numbers don't match
    if ~indext && isequal(sz2*r.ensembles(1),lsize) %%if ensemble missing
        o = o.*ones([sz,r.ensembles(1)]);    %%add ensemble
    else                                     %%numbers still don't match
        return                               %%return
    end                                      %%end if ensemble is missing
end                                          %%end if numbers don't match
o = reshape(o,[sz(1),lattice]);              %%reshape to lattice
end                                          %%end function