                                            
function e = Equilibrium()                          %%name of main function

in.name =      'Equilibrium spectrum';     %%name for simulation
in.points =    640;
in.ranges =    100;
in.noises =    [2,0];                      %%xnoises, knoises per point
in.ensembles = [100,5,1];                  %%samples,ensembles,parallel
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation 
in.da =        @(a,xi,~)  -a + xi(1,:)+1i*xi(2,:);   %%Derivative  
in.observe{1} =   @(a,~) a.*conj(a);       %%Observe  handle
in.observe{2} =   @(a,~) a.*conj(a);       %%Observe  handle
in.transforms = {0,1};
in.olabels =   {'|a(t)|^2','|a(w)|^2'};    %%labels for observables
in.compare =   {@(t,~) 1.+0*t, @(w,~)100./(pi*(1+w.^2))};%%Comparison handle

e            = xspde(in);                  %%Stochastic program
end                                             