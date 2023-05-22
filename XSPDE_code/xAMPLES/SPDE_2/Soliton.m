function [e]  =  Soliton()
%   [e]  =  SOLITON() tests xSPDE for a nonlinear Schroedinger equation
%   Tests a two-dimensional partial differential equation for:
%   (1) Inputting the x-grid into the initial condition
%   (2) Nonlinear derivative function in a PDE
%   (3) Setting up the linear response for the interaction picture
%   (4) Compare functions for a 2D simulation
%   (5) Graphing results in 1D and 2D, including transverse graphs
%   (6) Default values for transverse range and points in x
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name         =  'Soliton example';                    %%Graph title
p.dimensions   =  2;                                    %%Space-time dimn
p.initial      =  @(v,p)         sech(p.x);             %%Initialisation
p.deriv        =  @(a,~,~)       1i*a.*(conj(a).*a);    %%Nonlinear deriv
p.linear       =  @(p)           0.5*1i*(p.Dx.^2-1.0);  %%Linear terms
p.compare      =  @(p)           sech(p.x);             %%Comparison handle
p.transverse   =  1;                                    %%Transverse graph
p.diffplot     =  1;                                    %%Transverse graph
e              =  xspde(p);                             %%main program
end                                                     %%end of function
