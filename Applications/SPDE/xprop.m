function [ac,bv]  =  xprop(ac,p) 
%   [ac,bv]  = XPROP(ac,p) computes one step of linear propagation. 
%   Input:  field cells ac, parameter structure p.
%   Output: field cells ac after propagation
%   Expects time variable p.t set to start of propagation
%   Treats periodic boundaries and ODEs only
%   Uses p.propagator{c,1} to propagate, in k-space if necessary
%   For SDE cases, the propagator is a vector 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Needs: p.propagator, p.dimensions
%   Calls: fft,ifft
%   Called by: method
%   Licensed by P. D. Drummond, S. Kiesewetter (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bv = 0;                                          % initial boundary values
for c = 1:p.fieldcells                           % loop over field cells
  if ~isequal(p.propagator{c},1)                 % if propagator not unity
    pat = 0;                                     % set the patch to zero
    if p.setboundcell{c}                         % if boundary values
      [a,pat,sw,patlp] = xpatch(ac,c,p);         % get cell c patched field
      if sw == 0                                 % if sw = 0, use patlp
        pat = pat + patlp;                       % add patlp to patch
      else                                       % otherwise ignore patlp
        for d = 2:p.dimensions                   % loop over dimensions
          pat = pat + p.da{c,d}.*DS(pat,d,c,p);  % compute corrected patch
        end                                      % end loop over dimensions
      end                                        % end if sw = 0
    else                                         % else no boundary values
      a = ac{c};                                 % get cell c field
    end                                          % end if boundary values
    a = xtransform(a,c,0,p);                     % transformed input field
    a = p.propagator{c}.*a;                      % propagate in k-space
    ac{c} = xtransform(a,c,1,p)+pat;             % inverse transform field
  end                                            % end if propagator
  if p.quantum == 1                              % if Schrodinger wavefn
    ac{c} = ac{c}./sqrt(sum(conj(ac{c}).*ac{c},1:p.nfields));% normalise it
  elseif p.quantum == 2                          % if density matrix
    ac{c} = ac{c}./trace(ac{c});                 % project density matrix
  end                                            % end tests for quantum
end                                              % end loop over cells c
if p.setboundaries                               % if boundaries needed
    %p.t = p.t + p.dtr/(p.ipsteps);               % set end boundary time
    [ac,bv] = p.setbound(ac,p);                  % set boundary values
end                                              % end if boundaries
end                                              % end xprop function