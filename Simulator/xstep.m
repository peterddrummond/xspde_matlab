
function a = xstep(a0,a1,w,dt,p)                 % Default single step

%   in = XSTEP(a0,a1,w,dt,p) sets default integration step.
%   Input: initial a0, deriv fields a1, noise w, stepsize dt, parameters p
%   Output: cell array of integrated output fields
%   Called by: method
%   Needs: d.ca - dimensions of fields in cell array
%   Calls: deriv
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for c = 1:p.fieldcells                         % loop over field cells
    a{c} = a0{c} + reshape(p.deriv{c}(a1{:},w{:},p)*dt,p.d.ca{c});% derivs      
  end                                            % end field cells loop
end