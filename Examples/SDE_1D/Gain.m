function [e]  =  Gain()
%   e  =  GAIN() tests xSPDE for an SDE with sequential loss then gain
%   Tests a one-dimensional linear stochastic differential equation for:
%   (1) An initial SDE whose solution should be time-invariant
%   (2) Inputting two steps per point for accuracy in only the second case
%   (3) Using two noises in ordinary space
%   (4) Setting the second input equal to the first, with default transfer
%   (5) Changing the second structure by overwriting with new data
%   (6) Running xspde with an input cell array having two structures
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name       =  'Loss with noise';                 %%name for simulation
p.ranges     =  4;                                 %%ranges: t
p.noises     =  2;                                 %%xnoises per point
p.ensembles  =  [10000,1,10];                      %%ensembles
p.initial    =  @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation
p.deriv      =  @(a,w,r)   -a + w(1,:)+1i*w(2,:);  %%Derivative
p.observe    =  {@(a,~) a.*conj(a)};               %%Observe  handle
p.olabels    =  {'|a|^2'};                         %%labels for observables
p.compare    =  {@(p) 1+0*p.t};                    %%Comparison handle
p2           =  p;                                 %%Second input
p2.steps     =  2;                                 %%Steps/plotted point
p2.name      =  'Gain with noise';                 %%name for simulation
p2.deriv     =  @(a,w,~)  a + w(1,:)+1i*w(2,:);    %%Derivative
p2.compare   =  {@(p) 2*exp(2*(p.t-4))-1};         %%Comparison handle
e             =  xspde(p,p2);                    %%Stochastic program
end
