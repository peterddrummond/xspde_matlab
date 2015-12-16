function e = Kubo()                            %%XSPDE Kubo oscillator
cd ~
in.name =       'Kubo oscillator';             %%Name of simulation
in.ensembles =  [2000,1,4];                    %%samples,ensembles,parallel
in.steps =      2;                             %%Number of steps
in.step =       @xRK4;                         %%Use RK4 integrator
in.initial   =  @(w,r)     1+0*w ;             %%Initialisation  handle
in.da        =  @(a,z,r)  1i*z.*a ;            %%Derivative  handle
in.file  =      'Kubo.mat';                    %%Output filename
in.compare   =  {@(t,~,~) exp(-t/2)};          %%Comparison function
[e,in,data]       =  xsim(in);                 %%Stochastic program
in.name =      'Kubo2 oscillator'; 
e            =  xgraph(in);                    %%Graphics program
end
