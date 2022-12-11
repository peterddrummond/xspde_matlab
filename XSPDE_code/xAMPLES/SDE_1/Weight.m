function [e]  =  Weight()
%   [e]  =  WEIGHT() tests xSPDE for a weighted SDE
%   Tests for:
%   (1) weighted propagation
%   (2) breeding using breedeps
%   (3) difference plots
%   (4) second-order extrapolation
%   (5) using xcheck
%   (6) returning breed fraction with breedfrac
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name         =  'Weight';                        %%Name of simulation
p.ensembles    =  [10000,10,12];                   %%samples,ensembles
p.fields       =  2;                               %%Fields plus weights
p.points       =  9;                               %%Integration points
p.order        =  2;                               %%Extrapolation order
p.thresholdw   =  0.1;                             %%Breeding threshold
p.diffplot     =  1;                               %%Flag to plot diff.
p.initial      =  @(w,~) [1+w(1,:);0*w(2,:)];      %%Initialisation  handle
p.deriv        =  @(a,w,~)  [-a(1,:)+w(1,:);-a(2,:)+w(2,:)];%%Deriv
p.observe{1}   =  @(a,~) a(1,:);                   %%Observe  handle # 1
p.observe{2}   =  @(a,p) p.breedw;                 %%Observe  handle # 2
p.compare{1}   =  @(p) exp(-p.t);                  %%Comparison function
p.olabels{1}   =  '<a>';                           %%y label # 1
p.olabels{2}   =  '<fractional breeds per step>';  %%y label # 2
e              =  xcheck(2,p);                     %%Convergence checks
end
