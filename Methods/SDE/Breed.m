function [varargout] = Breed(varargin)
%   [a,nbreed] = Breed(a,p) gives a trajectory breeding step
%  Input weights a{p.fieldcells} are real or complex. The procedure is
%   - get the average absolute weight 
%   - remove any absolute weight below Wt = thresholdw*W
%   - copy maximum-weight trajectory at its place, divide both weights by 2
%   - repeat until all trajectories have relative weight >= threshold
%   Input:  fields 'a', parameters p with threshold 'p.thresholdw'. 
%   Output: modified fields 'a', fractional breeds: 'nbreed'.
%   Called by xpath
%   Licensed by Peter D. Drummond, (2023) - see License.txt, XSDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE NUMBERS, TAKE ABSOLUTE VALUE OF WEIGHTS  

nbreed  = 0;                                     %%breed counter
dobreed = 1;                                     %%breed switch
csize = nargin-1;                                %%input field size  
a = varargin(1:csize);                           %%input fields 
p = varargin{end};                               %%input parameters
wtindx = p.fieldcells;                           %%weight index
sz = cell(1,csize);                              %%initialize array size
szm = cell(1,csize);                             %%initialize matrix size
for c = 1:csize                                  %%loop over cells  
  sz{c} = size(a{c});                            %%array size
  szm{c}  = [prod(1:sz{c}(end-1)),sz{c}(end)];   %%matrix size
end                                              %%end loop over cells  
ne   = sz{wtindx}(end);                          %%get ensemble number

% BREED INDIVIDUAL TRAJECTORIES INSIDE EACH PARALLEL SET 

while (dobreed)                                  %%repeat while dobreed = 1
    w    = exp(real(a{wtindx}(1,:)));            %%absolute weights of set
    wm   = mean(w);                              %%mean absolute weight
    idx  = w/wm < p.thresholdw;                  %%weights below threshold
    if ~any(idx)                                 %%if none identified
      dobreed    = false;                        %%unset breed flag, exit
    else                                         %%if weights identified
      [~,maxi]   = max(w);                       %%identify max weight
      [~,mini]   = min(w);                       %%identify min weight 
      for c = 1:csize                            %%loop over cells
        if ~isequal(c,wtindx)                    %%if not weight index
          a{c}  = reshape(a{c},szm{c});          %%reshape to matrix
          a{c}(:,mini)  = a{c}(:,maxi);          %%replace min by max
          a{c}  = reshape(a{c},sz{c});           %%reshape to input
        end                                      %%end if not weight index
      end                                        %%end loop over cells
      a{wtindx}(1,mini) = a{wtindx}(1,maxi)/2;   %%set new weight to half
      a{wtindx}(1,maxi) = a{wtindx}(1,mini);     %%set old weight to half 
      nbreed = nbreed + 1;                       %%increment counter
    end                                          %%end if none identified
end                                              %%end repeat while dobreed

% RETURN FIELDS
                                          
nbreed  = nbreed/ne;                             %%rate of breeding
varargout = a;
varargout{nargin} = nbreed;
end                                              %%end xbreed function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XBREED FUNCTION