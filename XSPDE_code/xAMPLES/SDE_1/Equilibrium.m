function [e] = Equilibrium()
%   e  =  EQUILIBRIUM() tests xSPDE for an SDE with a known spectrum
%   Tests a one-dimensional linear stochastic differential equation for:
%   (1) Increasing the number of time points and range
%   (2) Changing the random number seed
%   (3) Setting the transform inputs to give a spectrum in frequency
%   (4) Setting two observables, one in time and one in frequency
%   (5) Plotting data in both time and frequency consecutively
%   (6) Setting comparisons using a cell array of inline functions
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'Equilibrium spectrum';             %%name for simulation
p.points       =  200;                                %%points in time
p.steps        =  2;                                  %%points in time
p.ranges       =  100;                                %%range in time
p.seed         =  241;                                %%set the random seed
p.noises       =  [2,0];                              %%xnoises per point
p.ensembles    =  [100,1,12];                         %%samples,ensembles
p.initial      =  @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2);  %%Initialisation
p.deriv        =  @(a,w,~)  -a + w(1,:)+1i*w(2,:);    %%Derivative
p.observe{1}   =  @(a,~) a.*conj(a);                  %%Observe  handle
p.observe{2}   =  @(a,~) a.*conj(a);                  %%Observe  handle
p.transforms   =  {0,1};                              %%Get frequency
p.olabels      =  {'|a(t)|^2','|a(\omega)|^2'};       %%label observables
p.compare      =  {@(p) 1., @(p)p.ranges(1)./(pi*(1+p.w.^2))};%%Compare
e              =  xspde(p);                           %%Stochastic program
end                                                   %%end of function
