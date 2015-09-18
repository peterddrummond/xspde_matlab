function [e,in] = Wiener()                      %%Wiener process
cd ~
 in.name =           'Wiener process';
 in.noises =          1;
 in.ensembles =      [1000,10];
 in.da =             @(a,z,r) z;               %%Derivative function
 in.raw =            1;                          %%Derivative function
 in.Matfile =        'Wiener';
 in.observe =        {@(a,r) a.^2};              %%Observable function
 in.olabels   =      {'<a^2>'};                  %%labels
 in.compare   =      {@(t,~) t};                 %%Comparison handle
 e        =          xspde(in);
end                                              %%end of main function
