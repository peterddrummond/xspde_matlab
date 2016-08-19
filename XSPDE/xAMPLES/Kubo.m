function e = Kubo()                            
%   e  =  KUBO() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Parallel ensemble operation
%   (2) Initial conditions for a local ensemble
%   (3) Complex multiplicative noise derivative using RK4 method
%   (4) Computing the simulation with xsim
%   (5) Changing the graph-name after the simulation
%   (6) Graphing data using a stored file-name in the 'in' structure
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
% cd ~                                         %%If needed to set directory
in.name =       'Kubo oscillator';             %%Name of simulation
in.ensembles =  [2000,4,1];                    %%samples,ensembles,parallel
in.fields = [1,0];
in.step =       @xRK4;                         %%Use RK4 integrator
%in.graphs =     1;
in.initial   =  @(w,r)     1+0*w ;             %%Initialisation  handle
in.da        =  @(a,z,r)  1i*z.*a(1,:) ;       %%Derivative  handle
in.file  =      'Kubo.mat';                    %%Output filename
in.compare   =  {@(t,~,~) exp(-t/2)};          %Comparison function
[e1,in,~]       =  xsim(in);                   %%Stochastic program
in.name =      'Kubo oscillator edited title';
e            =  xgraph(in.file,in);            %%Graphics program

end
