function [e] = Kubo()
%   [e]  =  KUBO() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Parallel ensemble operation
%   (2) Initial conditions for a local ensemble
%   (3) Complex multiplicative noise derivative using RK4 method
%   (4) Computing the simulation with xsim
%   (5) Changing the graph-name after the simulation
%   (6) Graphing data using a stored file-name in the 'in' structure
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name      =  'Kubo oscillator';              %%Name of simulation
p.ensembles =  [1000,8,1];                    %%samples,ensembles,parallel
p.fields    =  1;
p.method    =  @RK4;                           %%Use RK4 integrator
p.initial   =  @(w,p)     1+0*w ;              %%Initialisation  handle
p.deriv     =  @(a,w,p)  1i*w.*a(1,:) ;        %%Derivative  handle
p.file      =  'Kubo.mat';                     %%Output filename
p.compare   =  @(p) exp(-p.t/2);               %%Comparison function
e           =  xsim(p);                        %%Stochastic program
p.name      =  'Kubo oscillator edited title';
xgraph(p.file,p);                              %%Graphics program
end
