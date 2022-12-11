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

p.name         =  'Characteristic';                   %%simulation name 
p.dimensions   =  2;                                  %%dimension: 2 
p.initial      =  @(w,p)        sech(2.*(p.x+2.5));   %%Initialisation 
p.linear       =  @(p)         -p.Dx;                 %%First derivative
p.points       =  [51,35];                            %%lattice
p.olabels      =  {'a'; 'a^2'};                       %%observe labels
p.observe{1}   =  @(a,~)         a ;                  %%First observe
p.observe{2}   =  @(a,~)         a.^2;                %%next observe
p.axes{2}      =  {0,18};                             %%Axes for plots
p.compare{2}   =  @(p)        sech(2.*(p.t-2.5)).^2;  %%Compare handle
e              =  xspde(p);                           %%main program
end                                                   %%end of function
