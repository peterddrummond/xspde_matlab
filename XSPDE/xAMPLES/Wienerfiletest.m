function e = Wienerfiletest()                      
%   e  =  WIENER() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for: 
%   (1) Adding an input name for graphs
%   (2) Setting the noises to 1 and integration method to xEuler
%   (3) Setting ensembles(2) to get the sampling error
%   (4) Inputting an inline derivative
%   (5) Storing raw data in a file: read using 'load Wiener.mat' 
%   (6) Inputting observe functions, olabels and comparisons
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

cd ~                                             %%sets user directory
 in.name =           'Wiener process';           %%name of simulation
 in.step   =         @xEuler;                    %%first integration method
 in.noises =          1;                         %%noises set to one
 in.ensembles =      [1000,10];                  %%ensembles for averaging
 in.da =             @(a,z,r) z;                 %%Derivative function
 in.raw =            1;                          %%raw switch
 in.file =           'Wiener.mat';               %%name of file output
 in.observe =        {@(a,r) a.^2};              %%Observable function
 in.olabels   =      {'<a^2>'};                  %%labels
 in.compare   =      {@(t,~) t};                 %%Comparison handle
 [e,in]        =          xsim(in);              %%Runs xsim simulation
 e            =      e+xgraph(in.file);
end                                              %%end of main function
