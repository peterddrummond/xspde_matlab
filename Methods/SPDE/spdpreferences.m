 function p = spdpreferences(p)
%   p = SPDPREFERENCES(p) chooses default values for spde data.
%   Input:  input parameter structure,'p'
%   Output:  parameter structure, with default spde values set.
%   Called by: xpreferences
%   Needs:     xprefer
%   xSDE functions are licensed by Peter D. Drummond, (2024) - see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET SPDE PREFERENCES

nd =  p.dimensions;                              % number of dimensions
nc =  p.fieldcells;                              % number of field cells
p.numberaxis = xprefer(p,'numberaxis',1,0);      % switch for axis numbers 
p.boundfun =   xprefer(p,'boundfun',1,@xboundfun);
p.setbound =   xprefer(p,'setbound',1,@xsetbound);
p.boundaries = xprefer(p,'boundaries',1,{zeros(p.fields{1},2)});
p.dkp = pi./((p.points{1} - 1).* p.dx);          % DCT and DST step
p.dkv = prod(p.dk(2:nd));                        % k-space volume
p.kranges = p.dk.*p.points{1};                   % ranges in k-space
p.dkv = prod(p.dk(2:nd));                        % k-space volume
p.nspace = prod(p.points{1}(2:nd));              % Transverse lattice size
p.v =   p.dv*p.nspace;                           % lattice volume
p.kv =  p.dkv*p.nspace;                          % k-space volume
p.points2 = floor(p.points{1}/2);
p.plotpts = ((0:p.points{1}(1)-1)-floor((p.points{1}(1))/2));
p.boundaries{nc+1,nd} = [];
if ~isfield(p,"boundval")
    p.boundval = cell(nc,nd);
end
p.boundval{nc,nd+1} = [];
p.setboundaries = 0;                             % initialize boundary
for d = 2:nd
  pts = p.points{1}(d);
  p0  = 0:pts-1;                                 % index vector
  p1  = (1:pts-1);                               % truncated index vector
  p.xc{d} = p.origins(d) + p0*p.dx(d);           % n-th x-coords
  p.kc{d} = (p0-p.points2(d))*p.dk(d);           % n-th k-coords
  ac = cell(1,p.fieldcells);
  p.setboundcell = cell(1,p.fieldcells);
  for c = 1:p.fieldcells
    if isequal(p.boundval{c,d},[])
      ac{c} = zeros(p.d.ca{c});
      p.boundval{c,d} = xboundfun(ac,c,d,p);
    end  
    if isequal(p.boundaries{c,d},[])
        p.boundaries{c,d} = zeros([p.fields{c},2]);
    end   
    if ~isequal(p.boundaries{c,d},zeros([p.fields{c},2]))
      p.setboundaries = 1;                       %  Setboundaries flagged
      p.setboundcell{c} = 1;                     %  Setboundary(c) flagged
    else
      p.setboundcell{c} = 0;                     %  Setboundar(c) unflagged
    end
    b = p.boundaries{c,d};                       % store boundary switch 
    type  = b(:,1)+2*b(:,2);                     % get integer switch
    for j = 1:p.fields{1}                        % loop over fields
      switch type(j)                             % switch on boundary type
        case -3                                  % Robin boundaries
          p.kcp{c,j}{d} = p0*p.dkp(d);
        case -1                                  % Dirichlet-Robin
          p.kcp{c,j}{d} = [0,p1-1/2]*p.dkp(d);
        case 0                                   % periodic boundaries
          p.kcp{c,j}{d} = ifftshift(p.kc{d});    % n-th propagation k-coord
        case 1                                   % Robin-Dirichlet
          p.kcp{c,j}{d} = [p1-1/2,0]*p.dkp(d);
        case 3                                   % Dirichlet boundaries
          p.kcp{c,j}{d} = p0*p.dkp(d);
        otherwise                                % invalid boundary type
          error('Invalid boundary type: [%d,%d]',b(j,1),b(j,2));
      end                                        % end switch
    end
  end                                            % end loop over cells
end                                              % end loop over dimension
end                                              % end function 

function b = xboundfun(ac,c,d,p)                  % default boundary values
%   b = XBOUNDFUN(ac,c,d,p) calculates the default boundary values
%   Input is a cell of fields, ac, and a cell index c
%   Output are lower and upper boundary fields at dimension d
%   Boundaries are set in dimension d > 1 at upper and lower boundaries
%   Requires: p.t (current time), p.origin, p.indext
%   Called by xprop

sz = size(ac{c});
s1 = sz(1);
sz(1) = 1;
sz(d) = 1;
b=cell(s1,2);
if ~isfield(p,'boundval') || isequal(p.boundval{c,d},[]) % test boundval 
    for i = 1:s1
        b{i,1} = zeros(sz);                      % set to zeros as default
        b{i,2} = zeros(sz);
    end
else                                             % exists so use previous
    b = p.boundval{c,d};                         % previous boundary values
end 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END DEFAULT  FUNCTIONS