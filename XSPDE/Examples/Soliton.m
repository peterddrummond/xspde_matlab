function [e] = Soliton()                                

in.name =       'NLS soliton';
in.dimension =  2;                                      %%dimension: 1-4 
in.initial =    @(w,r)         sech(r.x);               %%Initialisation 
in.da =         @(a,~,r)       i*a.*(conj(a).*a);       %%Derivative 
in.linear =     @(D,r)         0.5*i*(D.x.^2-1.0);      %%laplacian
in.compare   =  @(t,~) 1+0*t;                           %%Comparison handle
[e,data] =      xspde(in);                           %%main program
end                                                     %%end of function
