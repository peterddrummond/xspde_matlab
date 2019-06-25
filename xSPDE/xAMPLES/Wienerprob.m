function e = Wienerprob()                      
%   e  =  WIENER() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for: 
%   (1) Storing raw data in a file: read using 'load Wiener.mat' 
%   (2) Setting the noises to 1 and integration method to xEuler
%   (3) Setting ensembles(2) to get the sampling error
%   (4) Inputting an inline derivative
%   (5) Using vector observe functions, olabels and comparisons
%   (6) Defining a second graph using the in.function transform
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 in.name =           'Wiener process distributions';%%name of simulation
 in.step   =         @xEuler;                    %%first integration method
 in.noises =          1;                         %%noises set to one
 in.ensembles =      [1000,10];                  %%ensembles for averaging
 in.initial =        @(z,r) z/2;                 %%Initial function
 in.da =             @(a,z,r) z;                 %%Derivative function
 in.observe{1} =     @(a,r) [a];  %%Observable function
 in.observe{2} =     @(a,r) [a;a.^2;a.^3;a.^4/30];  %%Observable function
 in.compare{2}   =   @(t,~) [0*t;.25+t;0*t;.1*(.25+t).^2];%%Comparisons
 in.transverse{1} =   5;                         %%Observable function
 in.olabels{1}   =    'a';                     %%labels
 in.probability{1} =   -5:.1:5;
 in.olabels{2}   =    '<a^n>';
 
 e        =          xspde(in);                  %%Runs xspde simulation
end                                              %%end of main function
