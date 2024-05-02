function [e]  =  SolitonDerivS()
%   e  =  SOLITONDERIVS() solves a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative
%   (2) Using additional points 
%   (3) Adding intermediate steps for higher accuracy
%   (4) Using periodic boundaries with finite differences
%   (5) Using derivatives to evaluate the knietic energy
%   (6) Using two sequential integrations with increased nonlinearity
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'SolitonDerivS';
p.dimensions   =  2;
p.order        =  4;
p.points       =  [101,51];
p.steps        =  10;
p.initial      =  @(v,p)   sech(p.x);                   %%Initialisation
p.deriv        =  @(a,w,p) 1i*a.*(conj(a).*a)+0.5*1i*(D2(a,p)-a);
p.observe{1}   =  @(a,p)   a.*conj(a);
p.observe{2}   =  @(a,p)   Int(abs(D1(a,p)).^2,p);
p.olabels      =  {'|a|^2','\int |da/dx|^2 dx'};
p1             =  p;
p1.deriv       =  @(a,w,p) 2.*1i*a.*(conj(a).*a)+0.5*1i*(D2(a,p)-a);
p1.name        =  'NLS using double nonlinearity';
e              =  xspde(p,p1);                        %%main program
end                                                     %%end of function
