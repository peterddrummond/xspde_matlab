function [e]  =  Weight()
%   [e]  =  WEIGHT() demonstrates a weighted SDE
%   Tests for:
%   (1) weighted propagation
%   (2) breeding by removal of low weights
%   (3) difference plots
%   (4) using multiple fields for weights
%   (5) using p.thresholdw
%   (6) returning breed fraction with breedw
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.name         =  'Weight';                      %%Name of simulation
p.ensembles    =  [10000,10,12];                  %%samples,ensembles
p.fields       =  {1,1};                         %%Fields plus weights
p.points       =  9;
p.order        =  2;
p.thresholdw   =  0.1;                           %%Breeding threshold
p.diffplot     =  1;                             %%Flag to plot diff.
p.initial      =  {@(w,~,~) 1+w,@(~,~,~) 0};     %%Initialisation  handle
p.deriv{1}     =  @(a,o,w,v,~) -a+w;             %%Derivative 1
p.deriv{2}     =  @(a,o,w,v,~) -o+v;             %%Derivative 2
p.expect{1}    =  @(a,o,~) a;                    %%Observe  handle # 1
p.observe{2}   =  @(~,~,p) p.breedw;             %%Observe  handle # 2
p.compare{1}   =  @(p) exp(-p.t);                %%Comparison function
p.olabels{1}   =  '<a>';                         %%y label # 1
p.olabels{2}   =  '<fractional breeds per step>';%%y label # 3
e              =  xcheck(2,p);                   %%
end