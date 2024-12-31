function w = xnoise(w1,nc,s,p)
%   w = XNOISE(w1,nc,s,p) returns p.noisecells + p.unoisecells noise arrays
%   Noise arrays have p.d.noises{j} dimensional noise fields in cell j.
%   Returns: p.noisecells of Gaussian noise * sqrt(s*p.noisefactor{j})
%            p.unoisecells of uniform noise on [0,1] * (s*p.noisefactor{j})
%
%   Where:   w1 is a previous noise term, used for error-checking
%            nc is a check index: nc = 2 for coarse steps, 
%            (only used if checks(1) = 1)
%            s is a time-step dependent scale factor for noise terms
%            p is the parameter structure
%            p.noisefactor{j} is a grid-dependent x-space scale factor
%         
%   If there are k-noises or uniform noises these are added as well.
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

w = cell(1,sum(p.noisecells));                   % initialize output cells
w{1} = 0;                                        % Initialize first to zero
j = 1;
for i = 1:p.noisecells(1)                        % Loop over noise cells
  if p.d.noise{i}(1)>0 
    w{j} = randn(p.d.noise{i});                  % Gaussian noises
    w{j} = w{j}*sqrt(s*p.noisefactor{i});        % Gaussian scale factor
    if ~isequal(w1,[]) && p.checks(1) && nc == 2 % If old noise & checks
      w{j} = (w{j}+w1{j})/2;                     % Use the average
    end                                          % End if old noise
  end
      j = j+1;
end                                              % End loop over gnoise
for i = 1:p.noisecells(2)                        % Loop over knoise cells
  if p.d.knoise{i}(1)>0 
    zk = randn(p.d.knoise{i});                   % Gaussian k-noise
    zk = p.nfilter{i}(zk,p);                     % knoise filter
    zk = xift(zk,p)*sqrt(s*p.knoisefactor{i});   % inverse FFT
    w{j} = zk;                                   % add to noise    
    if ~isequal(w1,[]) && p.checks(1) && nc == 2 % If old noise & checks
      w{j} = (w{j}+w1{j})/2;                     % Use the average
    end                                          % End if old noise
  end
  j = j+1;
end                                              % End knoise loop
for i = 1:p.noisecells(3)                        % Loop over uniform cells
  if p.d.unoise{i}(1)>0 
    w{j} = rand(p.d.unoise{i})*s*p.noisefactor{i};% rescale 
    if ~isequal(w1,[]) && p.checks(1) && nc == 2 % If old noise & checks
      w{j} = min(w{j},w1{j});                    % Use the minimum
    end                                          % End if old noise
  end
  j = j+1;
end                                              % End uniform noise loop
end                                              % End noise function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%