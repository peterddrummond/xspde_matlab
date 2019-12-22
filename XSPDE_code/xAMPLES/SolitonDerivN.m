function [e] = SolitonDerivN()
%   e  =  SOLITONDERIVN() tests xSPDE for a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative, in.da
%   (2) Using Neumann or zero derivative boundary conditions
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.dimension =     2;                                     
in.points =        [101,40];                                    
in.steps =         10;
in.initial =       @(w,r)   sech(r.x);                 %%Initialisation
in.observe{1} =    @(a,r)   a.*conj(a);
in.observe{2} =    @(a,r)   xint (abs(xd(a,r.Dx,r)).^2,r);
in.olabels    =    {'|a|^2','\int |da/dx|^2 dx'};
in.name =          'NLS soliton using finite differences + Neumann';
in.boundaries{2} =    [-1,-1];
in.transverse  =   {3};
in.da        =    @(a,~,r)    1i*a.*(conj(a).*a)+0.5*1i*(xd2(a,2,r)-a);
e =xspde(in);                                 %%main program
end                                                     %%end of function