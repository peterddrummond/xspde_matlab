function a  =  xshape(a,c,indext,p)                 
%   a  =  XSHAPE(a,c,indext,p) shapes a field to standard xSPDE dimensions. 
%   Input: field  'a', cell index 'c', flag 'indext', parameters 'p'.
%   Checks if indext = 1  is set for possible time-index
%   Can shape either xSPDE fields (indext = 0) or averages (indext = 1)
%   Broadcasts fields over ensembles if needed
%   Leaves other inputs with non-standard total components unchanged.
%   First dimension(s): field indices (indext = 0) or line index (indext = 1)
%   if p.indext is zero, the last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sz  = size(a);                               % initial size
if p.dimensions == 1                         % check if ODE or SDE
    p.d.space{c} = [];                       % space is initialized to null
end                                          % end check if ODE or SDE
if indext                                    % if time index is set
    sz1 = sz(1);                             % set line index
    sz2 = sz(2:end);                         % set remaining indices
    lattice = p.points{c};                   % lattice is space-time
else                                         % if time index is not set
    sz1 = sz(1:p.nfields);                   % set field indices
    sz2 = sz(p.nfields+1:end);               % set remaining indices
    lattice = [p.d.space{c},p.ensembles(1)]; % lattice is space-ensembles
end                                          % end if time index is set
if isequal(sz2,lattice)                      % if lattice matches indices
    return                                   % return with no changes
end                                          % end if lattice OK
sz2 = prod(sz2);                             % combine remaining components
lsize = prod(lattice);                       % get total lattice size
if ~isequal(sz2,lsize)                       % if sizes don't match
    if ~indext && isequal(sz2*p.ensembles(1),lsize) % if ensemble missing
        a = a.*ones([sz,p.ensembles(1)]);    % add ensemble
    else                                     % if sizes still don't match
        return                               % return with no changes
    end                                      % end if ensemble is missing
end                                          % end if numbers don't match
a = reshape(a,[sz1,lattice]);                % reshape to lattice
end                                          % end function