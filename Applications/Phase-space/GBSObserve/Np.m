function N1 = Np(a,p)
% C = NP(a,p) generates mean photon number per channel
% Uses an ordered phase-space representation
% p.phase = 1,2,3 for P-representation, Wigner, Q repreentation
% Input:  field amplitudes a, parameter structure, p.
% Output: photon numbers per mode.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sm  = (p.phase-1.0)/2;                           % ordering correction
if p.phase == 1 && size(a,1) == 2*p.modes        % if positive P
  N1 = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);  % normal-ordered numbers
else                                             % else not +P 
  N1 = a.*conj(a)-sm;                             % add ordering correction
end                                              % end if positive P
end                                              % End pn function