function [e] = SolitonAv()
%   e  =  SOLITONAV() tests xSPDE for a nonlinear Schroedinger equation
%   Tests a two-dimensional partial differential equation for:
%   (1) Putting print to 0 to minimize output printing
%   (2) Averaging over the spatial lattice
%   (3) Setting up the linear response for the interaction picture
%   (4) Compare functions for a 2D simulation
%   (5) Graphing results in 1D and 2D
%   (6) Default values for transverse range and points in x
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'SolitonAv: averages and integration';
p.dimensions   =  2;                               %%dimension: 1-4
p.observe{1}   =  @(a,p)      Ave(real(a),p.dx,p);
p.observe{2}   =  @(a,p)      Int(real(a),p.dx,p);
p.initial      =  @(v,p)      sech(p.x);           %%Initialisation
p.deriv        =  @(a,~,~)    1i*a.*(conj(a).*a);  %%Derivative
p.linear       =  @(p)      0.5*1i*(p.Dx.^2-1.0);  %%Laplacian
p.compare      =  {@(p) pi/10+0*p.t,@(p) pi+0*p.t};
p.olabels      =  {'<<a(x)>>','\int<a(x)>dx'};     %%Labels
p.pdimension   =  {1,1};
e              =  xspde(p);                        %%main program
end                                                %%end of function
