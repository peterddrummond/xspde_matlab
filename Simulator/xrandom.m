function w = xrandom(p)
%   w = XRANDOM(p) generates initial random number arrays
%   Generates p.randcells  gaussian random fields, up to length of array
%   Fields are returned sequentially, with x-correlated fields first.
%   Fields generated in k-space are transformed back to x-space
%   All fields are scaled by the standard random scaling factors
%   Any k-random or unifrom random inputs are added, note that:
%       a) noisefactor{j} is a grid-dependent x-space scale factor
%       b) knoisefactor{j} is a grid-dependent k-space scale factor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:  parameter structure p 
%   Output: random cell array w
%   Called by: p.xpath
%   Licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = cell(1,sum(p.randcells));                    % Create cell array
w{1} = 0;                                        % Initialize to zero
for n = 1:p.randcells(1)                         % Loop on gaussian cells
  w{n} = randn(p.d.rand{n});                     % Gaussian randoms
  w{n} = w{n}*sqrt(p.noisefactor{n});            % Scale Gaussian noise
end                                              % End loop on gaussians
for n = 1:p.randcells(2)                         % Loop on k cells
    kv = randn(p.d.krand{n});                    % Get k-randoms
    kv = p.rfilter{n}(kv,p);                     % Filter k-randoms
    kv = sqrt(p.knoisefactor{n})*xift(kv,p);     % Inverse FF transform 
    w{n+p.randcells(1)} = kv;                    % Combine randoms
end                                              % End loop on krandoms
for j = 1:p.randcells(3)                         % Loop over uniform cells
  n = j+p.randcells(1)+p.randcells(2);           % offset the index
  w{n} = rand(p.d.urand{j})*p.noisefactor{j};    % rescale 
end                                              % End uniform noise loop
end                                              % End xrandom function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%