function [e] = Soliton()
%   e  =  SOLITON() tests xSPDE for a nonlinear Schroedinger equation
%   Tests a two-dimensional partial differential equation for:
%   (1) Inputting the x-grid into the initial condition
%   (2) Nonlinear derivative function in a PDE
%   (3) Setting up the linear response for the interaction picture
%   (4) Compare functions for a 2D simulation
%   (5) Graphing results in 1D and 2D
%   (6) Default values for transverse range and points in x
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'NLS soliton';
in.dimension =  2;                                      %%dimension: 1-4
in.step =@xRK4;
%in.axes{1}   =  {0,-1};
in.initial =    @(w,r)         sech(r.x);               %%Initialisation 
in.da =         @(a,~,r)       1i*a.*(conj(a).*a);      %%Derivative 
in.linear =     @(r)           0.5*1i*(r.Dx.^2-1.0);    %%laplacian
in.compare   =  @(t,~) 1+0*t;                           %%Comparison handle
e =             xspde(in);                              %%main program
end                                                     %%end of function
