function e = Wiener()                      
%   e  =  WIENER() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for: 
%   (1) Storing raw data in a file: read using 'load Wiener.mat' 
%   (2) Setting the noises to 1 and integration method to xEuler
%   (3) Setting ensembles(2) to get the sampling error
%   (4) Inputting an inline derivative
%   (5) Using vector observe functions, olabels and comparisons
%   (6) Defining a second graph using the in.function transform
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 in.name =           'Wiener process';           %%name of simulation
 in.step   =         @xEuler;                    %%first integration method
 in.noises =          1;                         %%noises set to one
 in.ensembles =      [1000,10];                  %%ensembles for averaging
 in.initial =        @(z,r) 2*z;                 %%Initial function
 in.da =             @(a,z,r) z;                 %%Derivative function
 in.raw =            1;                          %%raw switch
 in.file =           'Wiener.mat';               %%name of file output
 in.observe =        @(a,r) [a;a.^2;a.^3;a.^4];  %%Observable function
 in.compare{1}   =   @(t,~) [0*t;4+t;0*t;3*(4+t).^2];%%Comparison handle
 in.function{2} =    @(o,r) o{1}(2,:);           %%Space-time function
 in.compare{2}   =   @(t,~) 4+t;                 %%Comparison handle
 in.function{3} =    @(o,r) o{1}(4,:);           %%Space-time function
 in.xfunctions{3} =  {@(t,r) (4+t).^2};          %%Space-time function
 in.compare{3}   =    @(tau,~) 3*tau;            %%Comparison handle
 in.olabels   =      {'<a^n>','<a^2>','<a^4>'};  %%labels
 in.glabels{1}   =    {'t'};                     %%labels
 in.glabels{2}   =    {'t'};                     %%labels
 in.glabels{3}   =    {'\tau = (4+t)^2'};        %%labels
 e        =          xspde(in);                  %%Runs xspde simulation
end                                              %%end of main function
