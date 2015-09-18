function e = SHO()         %%Simple harmonic oscillator with noise

 in.name =      'Simple harmonic oscillator'      %%Name 
 in.initial =   @(~,~) 1;                         %%Initialisation 
 in.da =        @(a,~,r) 1i*a;                    %%Derivative function   
 in.compare =   {@(t,in) cos(t)};                 %%Compare  handle
 e =xspde(in);
end                                               %%end of main function
