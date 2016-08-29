function e = Equilibrium()
%   e  =  EQUILIBRIUM() tests xSPDE for an SDE with a known spectrum
%   Tests a one-dimensional linear stochastic differential equation for:
%   (1) Increasing the number of time points and range
%   (2) Changing the random number seed
%   (3) Setting the transform inputs to give a spectrum in frequency
%   (4) Setting two observables, one in time and one in frequency
%   (5) Plotting data in both time and frequency consecutively
%   (6) Setting comparisons using a cell array of inline functions
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =      'Equilibrium spectrum';           %%name for simulation
in.points =    640;                              %%points in time
in.ranges =    100;                              %%range in time
in.seed =      240;                              %%set the random seed
in.noises =    [2,0];                            %%xnoises per point
in.ensembles = [100,5,1];                        %%samples,ensembles
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2);%%Initialisation 
in.da =        @(a,z,~)  -a + z(1,:)+1i*z(2,:);  %%Derivative  
in.observe{1} =   @(a,~) a.*conj(a);             %%Observe  handle
in.observe{2} =   @(a,~) a.*conj(a);             %%Observe  handle
in.transforms = {0,1};                           %%Transform to frequency
in.olabels =   {'|a(t)|^2','|a(\omega)|^2'};     %%labels for observables
in.compare =   {@(t,~) 1.+0*t, @(w,~)(100.16)./(pi*(1+w.^2))};%%Comparisons
e            = xspde(in);                        %%Stochastic program
end                                              %%end of function                                           