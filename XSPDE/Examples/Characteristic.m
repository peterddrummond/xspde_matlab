function [e] = Characteristic()                          %%name of function
%% Tests a simple characteristic solution of a first order PDE
in.name      =  'Characteristic'
in.dimension =  2;                                       %%dimension: 1-4 
in.initial =    @(w,r)      sech(2.*(r.x+2.5));          %%Initialisation 
in.da =         @(a,z,r)  0*a;                           %%Derivative 
in.linear =     @(D,r)         -D.x;                     %%First derivative
in.compare =    {@(t,in)        sech(2.*(t-2.5))};       %%Compare handle
e             = xspde(in);                            %%main program
end                                                      %%end of function
