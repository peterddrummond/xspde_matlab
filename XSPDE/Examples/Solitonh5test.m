function [e] = Solitonh5test()
%   e  =  SOLITON() tests xSPDE for a nonlinear Schroedinger equation
%   Tests a two-dimensional partial differential equation for:
%   (1) Inputting the x-grid into the initial condition
%   (2) Nonlinear derivative function in a PDE
%   (3) Setting up the linear response for the interaction picture
%   (4) Compare functions for a 2D simulation
%   (5) Graphing results in 1D and 2D
%   (6) Default values for transverse range and points in x
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
cd ~
in.name =       'NLS soliton';
in.dimension =  4;
in.numberaxis  = 1;
in.file      =  'soliton.h5';
in.initial =    @(w,r)         sech(r.x{1})%.*sech(r.x{2});               %%Initialisation 
in.da =         @(a,~,r)       1i*a.*(conj(a).*a);      %%Derivative 
in.linear =     @(D,r)         0.5*1i*(D{1}.^2-1.0);     %%laplacian
in.compare   =  @(t,~) 1+0*t;                           %%Comparison handle
in.images = 2;
[e,in] =             xsim(in);                              %%main program
e =             xgraph(in);                              %%main program
end                                                     %%end of function
