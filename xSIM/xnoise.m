function w = xnoise(w1,nc,s,p)
%   w = XNOISE(w1,s,p) generates a total of p.noisecells of noise arrays
%   Arrays have p.d.noises{j} dimensions of delta-correlated noise fields
%   Returns: (1) p.normalcells of normal noise * sqrt(s*p.noisefactor{j})
%            (2) p.uniformcells of uniform noise * (s*p.noisefactor{j})          
%   Here: s is a time-step dependent scale factor for delta-correlations
%   If there are k-noise inputs defined, these are added to (1), where:
%         a)noisefactor{j} is a grid-dependent x-space scale factor
%         b)knoisefactor{j} is a grid-dependent k-space scale factor.
%   If w1 present, returns an average (Gaussian) or minimum (uniform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:  previous noise cell array w1, scale s, parameter structure p 
%   Output: new noise cell array w
%   Called by: p.xpath
%   xSDE functions are licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = cell(1,p.noisecells+p.unoisecells);          %% initialize output cells
for j = 1:p.noisecells                           %% Loop over noise cells
  w{j} = randn(p.d.noises{j});                    %% Gaussian noises
  w{j} = w{j}*sqrt(s*p.noisefactor{j});          %% Gaussian scale factor
  if p.noises{j}(2) > 0                          %% if k-space noise
    zk = randn(p.d.knoises{j});                  %% get k-noise
    zk = p.nfilter{j}(zk,p);                     %% knoise filter
    zk = xift(zk,p)*sqrt(s*p.knoisefactor{j});   %% inverse FFT
    w{j}=[w{j};zk];                              %% add to noise
  end                                            %% end if k-space noise      
  if ~isequal(w1,[]) && p.checks(1) && nc == 2   %% If previous noise
    w{j} = (w{j}+w1{j})/2;                       %% Take the average
  end                                            %% End if previous noise
end                                              %% End Gaussian noise loop
for j1 = 1:p.unoisecells                         %% Loop over uniform cells
  j = j1+p.noisecells;                           %% offset the index
  w{j+p.noisecells } = rand(p.d.noises{j})*s*p.noisefactor{j};
  if ~isequal(w1,[]) && p.checks(1) && nc == 2   %% If previous noise
    w{j} = min(w{j},w1{j});                      %% Take the minimum
  end                                            %% End if previous noise
end                                              %% End uniform noise loop
end                                              %% End noise function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%