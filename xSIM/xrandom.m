function w = xrandom(p)
%   w = XRANDOM(p) generates initial random number arrays
%   Generates p.randcells  gaussian random fields, up to length of array
%   Fields are returned sequentially, with x-correlated fields first.
%   Fields generated in k-space are transformed back to x-space
%   All fields are scaled by the standard random scaling factors
%   If there are k-random inputs defined, these are appended, where:
%         a)noisefactor{j} is a grid-dependent x-space scale factor
%         b)knoisefactor{j} is a grid-dependent k-space scale factor.
%   For error-checking, xrandom 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:  parameter structure p 
%   Output: random cell array w
%   Called by: p.xpath
%   Licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = cell(1,p.randomcells);
for n = 1:p.randomcells                          %% Loop on gaussian cells
  %if p.dimensions > 1 
  w{n} = randn(p.d.randoms{n});                  %% Gaussian randoms
  w{n} = w{n}*sqrt(p.noisefactor{n});            %% Scale Gaussian noise
  if p.inrandoms{n}(2) > 0                           %% Test for k randoms
    kv = randn(p.d.krandoms{n});                 %% Get k-randoms
    kv = p.rfilter{n}(kv,p);                     %% Filter k-randoms
    kv = sqrt(p.knoisefactor{n})*xift(kv,p);     %% Inverse FF transform 
    w{n} = [w{n};kv];                            %% Combine randoms
  end                                            %% End if krandoms
end                                              %% End random loop
end                                              %% End xrandom function