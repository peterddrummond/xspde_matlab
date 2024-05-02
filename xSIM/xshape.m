function a  =  xshape(a,c,indext,p)                 
%   a  =  XSHAPE(a,c,indext,p) shapes a field to standard xSPDE dimensions. 
%   Input: field  'a', cell index 'c', flag 'indext', parameters 'p'.
%   Checks if indext is set for possible time-index
%   Can shape either xSPDE fields or averages
%   Broadcasts over an ensemble if needed
%   Leaves fields with non-standard total components unchanged.
%   First dimension is the field or line index 
%   if r.indext is zero, the last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sz  = size(a);                               %%initial size
if p.dimensions == 1
    p.d.space{c} = [];
end
if indext                                    %%if time index is set
    sz1 = sz(1);
    sz2 = sz(2:end);                         %%size of indices > 1
    lattice = p.points{c};                   %%lattice is space-time
else                                         %%if time index is not set
    sz1 = sz(1:p.nfields);
    sz2 = sz(p.nfields+1:end);
    lattice = [p.d.space{c},p.ensembles(1)];
end                                          %%end if time index is set
if isequal(sz2,lattice)                      %%if lattice OK
    return                                   %%return
end                                          %%end if lattice OK
sz2 = prod(sz2);                             %%store remaining components
lsize = prod(lattice);                       %%total lattice size
if ~isequal(sz2,lsize)                       %%if numbers don't match
    if ~indext && isequal(sz2*p.ensembles(1),lsize) %%if ensemble missing
        a = a.*ones([sz,p.ensembles(1)]);    %%add ensemble
    else                                     %%numbers still don't match
        return                               %%return
    end                                      %%end if ensemble is missing
end                                          %%end if numbers don't match
a = reshape(a,[sz1,lattice]);                %%reshape to lattice
end                                          %%end function