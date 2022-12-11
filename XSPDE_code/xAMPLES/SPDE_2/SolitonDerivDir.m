function [e]  =  SolitonDerivDir()
%   [e]  =  SOLITONDERIVDir() tests xSPDE for a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative, in.da
%   (2) Using Dirichlet, zero field boundary conditions
%   (3) Adds initial and additive noise terms
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.dimensions      =  2;
p.points          =  [101,40];
p.steps           =  10;
p.N               =  0.01;
p.initial         =  @(v,p)   sech(p.x)+v*p.N;               %%Initialisation
p.observe{1}      =  @(a,p)   a.*conj(a);
p.observe{2}      =  @(a,p)   Int (abs(D1(a,p)).^2,p);
p.olabels         =  {'|a|^2','\int |da/dx|^2 dx'};
p.name            =  'NLS soliton using finite differences + Dirichlet';
p.boundaries{2}   =  [1,1];
p.deriv           =  @(a,w,p) 1i*a.*(conj(a).*a)+0.5*1i*(D2(a,p)-a)+w*p.N;
e                 =  xspde(p);                               %%main program
end                                                          %%end of function
