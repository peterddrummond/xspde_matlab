function [e]       =  Wiener()
%   [e]  =  WIENER() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for:
%   (1) Storing raw data in a file: read using 'load Wiener.mat'
%   (2) Setting the noises to 1 and integration method to xEuler
%   (3) Setting ensembles(2) to get the sampling error
%   (4) Inputting an inline derivative
%   (5) Using vector observe functions, olabels and comparisons
%   (6) Defining a second graph using the in.function transform
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 p.name            =  'Wiener process';            %%name of simulation
 p.method          =  @Euler;                     %%first method
 p.noises          =  1;                           %%noises set to one
 p.ensembles       =  [250,4,10];                  %%ensembles used
 p.initial         =  @(v,~) 2*v;                  %%Initial function
 p.deriv           =  @(a,w,~) w;                  %%Derivative function
 p.rawdata         =  1;                           %%raw switch
 p.file            =  'Wiener.mat';                %%name of file output
 p.observe         =  @(a,~) [a;a.^2;a.^3;a.^4];   %%Observable function
 p.compare{1}      =  @(p) [0*p.t;4+p.t;0*p.t;3*(4+p.t).^2];%%Comparison
 p.output{2}       =  @(o,~) o{1}(2,:);            %%Space-time function
 p.compare{2}      =  @(p) 4+p.t;                  %%Comparison handle
 p.output{3}       =  @(o,~) o{1}(4,:)  ;          %%Space-time function
 p.xfunctions{3}   =  {@(t,p) (4+t).^2};           %%Space-time function
 p.compare{3}      =  @(p) 3*(4+p.t).^2;           %%Comparison handle
 p.olabels         =  {'< a^n >','< a^2 >','< a^4 >'};   %%labels
 p.glabels{3}      =  {'\tau = (4+t)^2'};          %%labels
 e                 =  xspde(p);                    %%Runs xspde simulation
end                                                %%end of main function
