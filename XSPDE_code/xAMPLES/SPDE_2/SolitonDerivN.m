function [e]  =  SolitonDerivN()
%   [e]  =  SOLITONDERIVN() tests xSPDE for a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using a trigonometric propagator
%   (2) Using Neumann or zero derivative boundary conditions
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.dimensions      =  2;                                     
p.points          =  [101,101];
p.ranges          =  [10,15];
p.initial         =  @(v,p)   sech(p.x);                 %%Initialisation
p.observe{1}      =  @(a,p)   a.*conj(a);
p.observe{2}      =  @(a,p)   Int(abs(D1(a,2,p)).^2,p);
p.olabels         =  {'|a|^2','\int |da/dx|^2 dx'};
p.name            =  'SolitonDerivN: uses spectral methods + Neumann';
p.boundaries{2}   =  [-1,-1];
p.transverse      =  {3};
p.deriv           =  @(a,~,p)    1i*a.*(conj(a).*a);
p.linear          =  @(p) 0.5*1i*(p.Dx.^2-1);
e                 =  xspde(p);                           %%main program
end                                                      %%end of function
