function e = Kubo()                            %%XSPDE Kubo oscillator
cd ~
in.name =       'Kubo oscillator';
in.ensembles =  [2000,1,4];                  %%samples,ensembles,parallel
in.steps =      2;
in.step =       @xRK4;
in.initial   =  @(w,r)     1+0*w ;              %%Initialisation  handle
in.da        =  @(a,z,r)  1i*z.*a ;           %%Derivative  handle
in.file  =      'Kubo.h5';                   %%Comparison numbers?
in.compare   =  {@(t,~,~) exp(-t/2)};          %%Comparison function
e            =  xsim(in)                  %%Stochastic program
e            =  xgraph('',in)                 %%Stochastic program
end
