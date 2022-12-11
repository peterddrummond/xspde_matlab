function [a,breedw] = xbreed(a,r)
%   [a,breedfrac] = xbreed(a,r) carries out a trajectory breeding step
%  Each input matrix a is sampled with r.ensembles(1) vector samples
%  Input weights a(nb,:) are logarithmic, or even complex. The procedure is
%   - determine the average absolute weight over a(nb,:), where nb=r.fields
%   - determine any absolute weight below W=r.thresholdw*meanweight
%   - if so, remove the minimum absolute-weight trajectory
%   - copy maximum-weight trajectory at its place, divide weights by 2
%   - repeat until all trajectories have relative weight >= thresholdw
%   Note: this version is designed for weighted SDEs, not PSDEs
%   Input:  field 'a',  lattice 'r'.
%   Output: log weights in field 'a', fractional breeds per step: breedw.
%   Licensed by Simon Kiesewetter and Peter D. Drummond (2021), see License

% INITIALIZE NUMBERS, TAKE EXPONENTIAL OF WEIGHTS  

nbreed  = 0;                                     %%breed counter 
dobreed = 1;                                     %%switch for breeding
nw      = r.fields;                              %%set weight id

% BREED INDIVIDUAL TRAJECTORIES INSIDE EACH PARALLEL SET 

while (dobreed)                                  %%repeat while dobreed = 1
    w    = exp(real(a(nw,:)));                   %%absolute weights of set
    wm   = mean(w);                              %%mean absolute weight
    idx  = w/wm < r.thresholdw;                  %%weights below threshold
    if ~any(idx)                                 %%if none identified
      dobreed    = false;                        %%unset breed flag, exit
    else                                         %%if weights identified
      [~,maxi]   = max(w);                       %%identify max weight
      [~,mini]   = min(w);                       %%identify min weight 
      a(:,mini)  = a(:,maxi);                    %%replace min by max
      a(nw,mini) = a(nw,maxi)-log(2);            %%set new weight to half 
      a(nw,maxi) = a(nw,mini);                   %%set old weight to half 
      nbreed     = nbreed + 1;                   %%increment counter
    end                                          %%end if none identified
end                                              %%end repeat while dobreed

% RETURN LOGARITHM OF WEIGHTS  
                                          
a(nw,:) = log(a(nw,:));                          %%turn weights to logs
breedw  = nbreed/r.ensembles(1);                 %%fraction of breed events
end                                              %%end xbreed function