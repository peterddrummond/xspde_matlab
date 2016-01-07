                                            
function e = Equilibriumtest()                          %%name of main function

in.name =      'Equilibrium spectrum';     %%name for simulation
in.points =    3200;
in.ranges =    500;
in.noises =    [2,0];                      %%xnoises, knoises per point
in.ensembles = [1000,1,32];                  %%samples,ensembles,parallel
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation 
in.da =        @(a,z,~)  -a + z(1,:)+1i*z(2,:);   %%Derivative  
in.observe{1} =   @(a,~) a.*conj(a);       %%Observe  handle
in.observe{2} =   @(a,~) a.*conj(a);       %%Observe  handle
in.transforms = {0,1};
in.olabels =   {'|a(t)|^2','|a(w)|^2'};    %%labels for observables
in.compare =   {@(t,~) 1.+0*t, @(w,~)in.ranges(1)./(pi*(1+w.^2))};%%Comparison handle

e            = xspde(in);                  %%Stochastic program
end                                             