function [ac,varargout]  =  xsetbound(ac,p) 
%   ac = XSETBOUND(ac,p) sets boundary values for propagation. 
%   Input:  field cells ac, parameter structure p.
%   Output: field cells ac, after boundary setting, with boundary values
%   Imposes Dirichlet boundaries in space if needed.
%   This function requires p.nfields = 1: fields are scalars or vectors
%   Boundary types and values are set for all field indices,
%   to all transverse dimension boundaries, and can change dynamically.
%   Boundary types are defined through p.boundaries{c,d}(i,b) for cell c,
%   space dimension (d=2,3..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d>1 is the transverse space dimension index, and options are:
%   (a) p.boundaries{c,d}(i,b)  = 0 gives the default, periodic 
%   (b) p.boundaries{c,d}(i,b)  = -1 gives Neumann, prescribed derivative
%   (c) p.boundaries{c,d}(i,b)  = 1 gives Dirichlet, prescribed field
%   Boundary values are set on fields through boundval{c,d}, which is
%   dynamically returned by the function p.boundfun(a,c,d,p)
%   The boundval shape matches the field, except in the boundary direction
%   In the direction of the derivative, index d, one has b = 1
%   Values are set sequentially by the dimension (d).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Needs:     p.fieldcells, p.propagator, p.dimensions, p.boundaries,
%              p.setboundaries
%   Calls:     boundfun
%   Called by: xprop
%   Licensed by P. D. Drummond & S. Kiesewetter (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c = 1:p.fieldcells
  a = ac{c};                                     % get the c field 
  sz = size(a);                                  % get the input size
  space = prod(sz(2:p.dimensions));              % space size
  bval = cell(p.fieldcells,p.dimensions); 
  if  space > 1                                  % If SPDE  
    for d = 2:p.dimensions                       % loop on space dimensions
      bound = p.boundaries{c,d};                 % store boundary switch
      if ~isequal(bound,zeros(size(bound)))      % check for boundaries
        bval{c,d} = p.boundfun(ac,c,d,p);        % get the boundary values
        if any(bound(:) == 1)                    % if any Dirichlet
          bv = bval{c,d};                        % initialize boundaries
          for i=1:sz(1)                          % loop over field index
            ind = p.ind;                         % first index = row
            ind{1} = i;
            if bound(i,1) == 1                   % If lower Dirichlet
              ind{d} = 1;                        % index d = 1        
              a(ind{:}) = bv{i,1};               % Set lower boundary val.
            end                                  % end if lower
            if bound(i,2) == 1                   % If upper Dirichlet
              ind{d} = sz(d);                    % index d = last             
              a(ind{:}) = bv{i,2};               % Set upper boundary val.
            end                                  % end if upper
          end                                    % end loop on field index
        end                                      % end check if Dirichlet
      end                                        % end check if boundaries
    end                                          % end loop on dimension
    ac{c} = reshape(a,sz);                       % reshape to input size
  end                                            % end if SPDE
end                                              % end cell loop 
if nargout == 2
  varargout{1} = bval;
end
end                                              %% end function