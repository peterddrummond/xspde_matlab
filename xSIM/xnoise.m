function w = xnoise(w1,nc,s,p)
%   w = XNOISE(w1,nc,s,p) returns p.noisecells + p.unoisecells noise arrays
%   Noise arrays have p.d.noises{j} dimensional noise fields in cell j.
%   Returns: p.noisecells of Gaussian noise * sqrt(s*p.noisefactor{j})
%            p.unoisecells of uniform noise on [0,1] * (s*p.noisefactor{j})
%
%   Where:   w1 is a previous noise term, used for error-checking
%            nc is a check index: nc = 2 for coarse steps, if checks(1) = 1
%            s is a time-step dependent scale factor for noise terms
%            p is the parameter structure
%            p.noisefactor{j} is a grid-dependent x-space scale factor
%         
%   If there are k-noises (ie, p.noises{j}(2) > 0), these are added to (1),
%
%   Where:   p.knoisefactor{j} is a grid-dependent k-space scale factor.
%
%   If w1 is present, returns an average (Gaussian) or minimum (uniform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:  previous noise array w1, check index nc, scale s, parameters p 
%   Output: new noise cell array w
%   Called by: xpath
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = cell(1,p.noisecells+p.unoisecells);          % initialize output cells
w{1} = 0;                                        % Initialize to zero
for j = 1:p.noisecells                           % Loop over noise cells
  w{j} = randn(p.d.noises{j});                   % Gaussian noises
  w{j} = w{j}*sqrt(s*p.noisefactor{j});          % Gaussian scale factor
  if p.noises{j}(2) > 0                          % if k-space noise
    zk = randn(p.d.knoises{j});                  % Gaussian k-noise
    zk = p.nfilter{j}(zk,p);                     % knoise filter
    zk = xift(zk,p)*sqrt(s*p.knoisefactor{j});   % inverse FFT
    w{j}=[w{j};zk];                              % add to noise
  end                                            % end if k-space noise      
  if ~isequal(w1,[]) && p.checks(1) && nc == 2   % If old noise & checks
    w{j} = (w{j}+w1{j})/2;                       % Use the average
  end                                            % End if old noise
end                                              % End Gaussian noise loop
for j1 = 1:p.unoisecells                         % Loop over uniform cells
  j = j1+p.noisecells;                           % offset the index
  w{j+p.noisecells } = rand(p.d.noises{j})*s*p.noisefactor{j}; % rescale 
  if ~isequal(w1,[]) && p.checks(1) && nc == 2   % If old noise & checks
    w{j} = min(w{j},w1{j});                      % Use the minimum
  end                                            % End if old noise
end                                              % End uniform noise loop
end                                              % End noise function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%