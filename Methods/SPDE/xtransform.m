function a  =  xtransform(a,c,inv,p) 
%   a = XTRANSFORM(a,c,inv,p) returns spectral transforms for cell c
%   Output: a, the spectral transform for the selected bcs
%   The field is assumed periodic unless there are specified boundaries
%   Use inv = 0 for forward transform, inv = 1 for inverse
%   Transform types are defined through p.boundaries{c,d}(i,b) for cell c,
%   dimension d = 2,3.., field index i = 1,2.. and boundary b = 1,2 
%   Here d > 1 is the space dimension index, and options are:
%
%   (a) p.boundaries{c,d}(i,b)  = 0  gives the default, periodic 
%   (b) p.boundaries{c,d}(i,b)  = -1 gives Neumann, prescribed derivative
%   (c) p.boundaries{c,d}(i,b)  = 1  gives Dirichlet, prescribed field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: xprop
%   Needs: p.boundaries, p.ind
%   Calls: xdct1,xdct2,xdct3,xdst1,xdst2,xdst3,fft,ifft
%   Licensed by P. D. Drummond, S. Kiesewetter (2024): see License 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1 = p.ind;                                      % initialize index cells
sz = size(a);                                    % field size
sz1 = [1,sz(2:end)];                             % one component size
s1 = sz(1);                                      % internal index size
irange = 1:s1;                                   % initialize index range
type = zeros(s1,p.dimensions);                   % initialize type
for d = 2:p.dimensions                           % loop over dimensions
  if sz(d) > 1                                   % check size > 1   
    bnd  = p.boundaries{c,d};                    % store boundary type
    type(:,d) = bnd(irange,1)+2*bnd(irange,2);   % boundary type index
  end                                            % end check size > 1   
end                                              % end loop over dimensions
for i = irange                                   % loop over components
  i1{1} = i;                                     % initialize index cells
  at = reshape(a(i,:),sz1);                      % store field component
  for d = 2:p.dimensions                         % loop over dimensions
    if sz(d) > 1                                 % check size > 1   
      if inv == 0                                % check if forward                        
        switch type(i,d)                         % switch on boundary case
        case -3                                  % Robin-Robin boundaries
          at = xdct1(at,d);                      % cosine-1 transform
        case -1                                  % Dirichlet-Robin boundary     
          at = xdst3(at,d);                      % take sine-3 transform
        case 0                                   % periodic boundaries
          at = fft(at,[],d);                     % Fourier transform
        case 1                                   % Robin-atirichlet
          at = xdct3(at,d);                      % take cosine-3 transform
        case 3                                   % Dirichlet-Dirichlet         
          at = xdst1(at,d);                      % take sine-1 transform            
        otherwise                                % invalid boundary type
          error('Invalid boundary: [%d,%d]',bnd(i,1),bnd(i,2));
        end                                      % end switch     
      else                                       % check if inverse    
        switch type(i,d)                         % inverse type 
        case -3                                  % Robin-Robin boundaries
          at = xdct1(at,d);                      % cosine-1 transform
        case -1                                  % Dirichlet-Robin
          at = xdst2(at,d);                      % sine-2 transform
        case 0                                   % atirichlet-Robin
          at = ifft(at,[],d);                    % Fourier transform
        case 1                                   % Robin-atirichlet 
          at = xdct2(at,d);                      % cosine-2 transform
        case 3                                   % Dirichlet-Dirichlet       
          at = xdst1(at,d);                      % sine-1 transform
        end                                      % end switch
      end                                        % end if inverse
    end                                          % end if size
  end                                            % end dimensions loop
  a(i1{:}) = at;                                 % store transform
end                                              % end components loop
end                                              % end xtransform function                                            