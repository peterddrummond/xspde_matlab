function [e] = SolitonDerivS()
%   e  =  SOLITONDERIVS() solves a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative, in.da
%   (2) Using additional points and steps for higher accuracy
%   (3) Changing the integration method to midpoint with no linear term
%   (4) Using periodic boundaries with spectral derivatives
%   (5) Using derivatives to evaluate the knietic energy
%   (6) Using two sequential integrations
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =         'NLS soliton using spectral derivatives';
in.dimension =     2;                                     
in.points =        [101,51];                                    
in.steps =         10;
in.initial =       @(w,r)   sech(r.x);                 %%Initialisation
in.da =           @(a,w,r) 1i*a.*(conj(a).*a)+0.5*1i*(xd(a,r.Dx.^2,r)-a);
in.observe{1} =    @(a,r)   a.*conj(a);
in.observe{2} =    @(a,r)   xint (abs(xd(a,r.Dx,r)).^2,r);
in.olabels    =    {'|a|^2','\int |da/dx|^2 dx'};
in1    =           in;
in1.da =           @(a,w,r) 1i*a.*(conj(a).*a)+0.5*1i*(xd2(a,r)-a);
in1.name =         'NLS soliton using periodic finite differences';
e =xspde({in,in1});                                     %%main program
end                                                     %%end of function