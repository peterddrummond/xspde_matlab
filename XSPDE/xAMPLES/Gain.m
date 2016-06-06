function [e] = Gain()                      
%   e  =  GAIN() tests xSPDE for an SDE with sequential loss then gain
%   Tests a one-dimensional linear stochastic differential equation for:
%   (1) An initial SDE whose solution should be time-invariant
%   (2) Inputting two steps per point for accuracy in only the second case
%   (3) Using two noises in ordinary space
%   (4) Setting the second input equal to the first, with default transfer
%   (5) Changing the second structure by overwriting with new data
%   (6) Running xspde with an input cell array having two structures
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =      'Loss with noise';          %%name for simulation
in.ranges =    4;                          %%ranges: t
in.noises =    [2,0];                      %%xnoises, knoises per point
in.ensembles = [10000,8,1];                %%samples,ensembles,parallel
in.ensembles = [10000,1];                  %%samples,ensembles,parallel
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation 
in.da =        @(a,xi,r)  -a + xi(1,:)+1i*xi(2,:);   %%Derivative  
in.observe =   {@(a,~) a.*conj(a)};        %%Observe  handle
in.olabels =   {'|a|^2'};                  %%labels for observables
in.compare =   {@(t,~) 1+0*t};             %%Comparison handle
in2        =   in;                         %%Second input
in2.steps  =    2;                         %%Steps per plotted point
in2.name =     'Gain with noise';          %%name for simulation
in2.da =        @(a,xi,~)  a + xi(1,:)+1i*xi(2,:);   %%Derivative 
in2.compare =   {@(t,~) 2*exp(2*t)-1};     %%Comparison handle
e           =   xspde({in,in2});           %%Stochastic program
end                                             
