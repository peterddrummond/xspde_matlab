function [e] = Characteristic()                          %%name of function
%   e  =  CHARACTERISTIC() tests the exact solution of a first order PDE
%   Tests a two-dimensional linear partial differential equation for:
%   (1) An initial condition as a sech funtion with an offset
%   (2) Inputting no derivative term to use the default value
%   (3) An exact compare function with a nontrivial time dependence
%   (4) Demonstrates periodic boundary conditions
%   (5) Using two distinct compare functions
%   (6) Using a cell array of observe labels
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name      =  'Characteristic';                      %%simulation name 
in.dimension =  2;                                     %%dimension: 2 
in.initial =    @(w,r)        sech(2.*(r.x+2.5));      %%Initialisation 
in.print =2; 
in.linear =     @(r)         -r.Dx;                   %%First derivative
in.olabels =    {'a'; 'a^2'};                          %%observe labels
in.observe{1} = @(a,r)         a ;                      %%First observe
in.observe{2} = @(a,r)         a.^2;                   %%First observe
in.compare{1} = @(t,in)        sech(2.*(t-2.5));       %%Compare handle
in.compare{2} = @(t,in)        sech(2.*(t-2.5)).^2;    %%Compare handle
e             = xspde(in);                             %%main program
end                                                    %%end of function
