function [e] = Solitontest()                                

in.name =       'NLS soliton with numerical axis label';
in.dimension =  2;                                      %%dimension: 1-4 
in.numberaxis = 1;
in.initial =    @(w,r)         sech(r.x{1});               %%Initialisation 
in.da =         @(a,~,r)       i*a.*(conj(a).*a);       %%Derivative 
in.linear =     @(D,r)         0.5*i*(D{1}.^2-1.0);      %%laplacian
in.compare   =  @(t,~) 1+0*t;                           %%Comparison handle
[e,data] =      xspde(in);                              

%%main program
end                                                     %%end of function
