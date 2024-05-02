function [ac,bv]  =  xprop(ac,p) 
%   [ac,bv]  = XPROP(ac,p) computes one step of linear propagation. 
%   Input:  field cells ac, parameter structure p.
%   Output: field cells ac after propagation, boundary value cells bv
%   Treats combinations of Dirichlet, Robin/Neuman or periodic boundaries
%   Each dimension is individually set at either end of every interval
%   All field components are treated indivdually
%   Uses p.propagator{c} to propagate, in k-space if necessary
%   For SDE cases, the propagator is a vector 
%   If p.propagator{1} = 1, no action apart from setting boundaries
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Needs: p.setboundaries, p.propagator, p.dimensions, p.boundaries
%   Calls: p.setbound,xdct1,xdct2,xdct3,xdst1,xdst2,xdst3,fft,ifft
%   Called by: method
%   Licensed by P. D. Drummond & S. Kiesewetter (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bv = cell(p.fieldcells,p.dimensions);            % boundary value cell 
bc = cell(1,p.dimensions);                       % boundary correction cell
for c = 1:p.fieldcells                           % loop over field cells c  
  if ~isequal(p.propagator{c},1)                 % If propagator needed
    a = ac{c};                                   % store field
    sz = size(a);                                % field size
    s1  = sz(1);                                 % internal index size
    type = zeros(s1,p.dimensions);               % initialize type
    for d = 2:p.dimensions                       % loop over dimension
      b = num2cell(zeros(s1,2));                 % initial boundary value
      if sz(d) > 1                               % If (S)PDE in dimension d
        bnd  = p.boundaries{c,d};                % store boundary type
        type(:,d) = bnd(:,1)+2*bnd(:,2);         % boundary type index
        if any(type(:,d))                        % If an index not periodic
          b  = p.boundfun(ac,c,d,p);             % get boundary values                      
          xa = p.origins(d); dx = p.ranges(d);   % get origin and range
          x  = p.r{d} - xa;                      % relative x coordinate
        end                                      % end if not periodic      
        for i = 1:s1                             % loop over field index
          in = p.ind;                            % initialize index cells
          in{1} = i;                             % set first index to i
          b1 = 0;                                % initial correction term
          et = 0;                                % initial dynamic bdary
          switch type(i,d)                       % switch on boundary case
            case -3                              % Robin-Robin boundaries
              db = (b{i,2} - b{i,1})/dx;         % difference in boundaries
              b1 = b{i,1}.*x + 0.5*x.^2.*db;     % static bdary correction            
              et = p.da{c}(i,d)*db;              % dynamic bdary correction            
            case -1                              % Dirichlet-Robin boundary
              b1 = b{i,1} + x.*b{i,2};           % static bdary correction     
            case 1                               % Robin-Dirichlet
              b1 = b{i,2} + (x - dx).*b{i,1};    % static bdary correction    
            case 3                               % Dirichlet-Dirichlet
              db = (b{i,2} - b{i,1})/dx;         % difference in boundaries          
              b1 = b{i,1} + x.*db;               % Dirichlet-Dirichlet
          end                                    % end switch     
          a(in{:}) = a(in{:}) - b1;              % subtract correction  
          bc{i,d}  = b1 + et;                    % store correction 
        end                                      % end loop over fields
      end                                        % end if (S)PDE
      bv{c,d} = b;                               % store boundary value
    end                                          % end loop on dimension
    for d = 2:p.dimensions                       % loop over dimension
      if sz(d) > 1                               % If (S)PDE in dimension d
        for i = 1:s1                             % loop over field index
          in{1} = i;                             % set first index to i
          switch type(i,d)                       % switch on boundary case
            case -3                              % Robin-Robin boundaries              
              a(in{:}) = xdct1(a(in{:}),d);      % cosine-1 transform
            case -1                              % Dirichlet-Robin boundary
              a(in{:}) = xdst3(a(in{:}),d);      % take sine-3 transform
            case 0                               % periodic boundaries
              a(in{:}) = fft(a(in{:}),[],d);     % Fourier transform
            case 1                               % Robin-Dirichlet              
              a(in{:}) = xdct3(a(in{:}),d);      % take cosine-3 transform
            case 3                               % Dirichlet-Dirichlet
              a(in{:}) = xdst1(a(in{:}),d);      % take sine-1 transform            
            otherwise                            % invalid boundary type
              error('Invalid boundary: [%d,%d]',bnd(i,1),bnd(i,2));
          end                                    % end switch
        end                                      % end loop over fields
      end                                        % end if (S)PDE
    end                                          % end loop on dimension
    a = p.propagator{c}.*a;                      % propagate in k-space
    for d = 2:p.dimensions                       % loop over dimension
      if sz(d) > 1                               % If (S)PDE in dimension d
        for i = 1:s1                             % loop over field index
          in{1} = i;                             % set first index to i
          if type(i,d) == 0                      % periodic boundaries
            a(in{:}) = ifft(a(in{:}),[],d);      % Fourier transform
          else                                   % not periodic 
            switch type(i,d)                     % switch on boundary type
            case -3                              % Robin-Robin boundaries
              a(in{:}) = xdct1(a(in{:}),d);      % cosine-1 transform
            case -1                              % Dirichlet-Robin
              a(in{:}) = xdst2(a(in{:}),d);      % sine-2 transform
            case 1                               % Robin-Dirichlet 
              a(in{:}) = xdct2(a(in{:}),d);      % cosine-2 transform
            case 3                               % Dirichlet-Dirichlet       
              a(in{:}) = xdst1(a(in{:}),d);      % sine-1 transform
            end                                  % end switch
            a(in{:}) = a(in{:}) + bc{i,d};     % store propagated slice
          end                                    % end if periodic
        end                                      % end components loop
      end                                        % End if (S)PDE case
    end                                          % end loop on dimension
    if all(sz(2:p.dimensions) == 1)              % pure ODE case
      a = p.propagator{c}.*a;                    % linear propagate in time
    end                                          % end pure ODE case
    ac{c} = a;                       % restore field and size    
  end                                            % end if propagator
  if p.quantum == 1                              % if Schrodinger wavefn
    ac{c} = ac{c}./sqrt(sum(conj(ac{c}).*ac{c},1:p.nfields));% normalise it
  elseif p.quantum == 2                          % if density matrix
    ac{c} = ac{c}./trace(ac{c});                 % project density matrix
  end                                            % end tests for quantum
end                                              % end loop over cells c
if p.setboundaries                               % else if boundaries
  [ac,bv] = p.setbound(ac,p);                    % set boundary values
end
end                                              % end xprop function