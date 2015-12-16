function e = Nonlinear()                            %%XSPDE Kubo oscillator
in.name =       'Nonlinear oscillator, RK4';
in.ensembles =  [2000,1,6];                  %%samples,ensembles,parallel
in.steps =      2;
in.points =     100;
in.step       = @xRK4;
in.initial   =  @(w,r)   0.5+0.25*w;        %%Initialisation  handle
in.da        =  @(a,xi,r) (a-a.^3)+xi;         %%Derivative  handle
in.observe    =  {@(a,r) a.^2, @(a,r) a.^4, @(a,r) a.^6};
in.olabels    =  {'<a^2>','<a^4>','<a^6>'};
e            =  xspde(in);                  %%Stochastic program
end
