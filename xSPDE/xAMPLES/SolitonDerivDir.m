function [e] = SolitonDerivDir()
%   e  =  SOLITONDERIVDir() tests xSPDE for a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative, in.da
%   (2) Using Dirichlet or zero field boundary conditions
%   (3) Adds initial and additive noise terms
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.dimension =     2;                                     
in.points =        [101,40];                                    
in.steps =         10;
in.N =             0.01;
in.initial =       @(w,r)   sech(r.x)+w*r.N;               %%Initialisation
in.observe{1} =    @(a,r)   a.*conj(a);
in.observe{2} =    @(a,r)   xint (abs(xd(a,r.Dx,r)).^2,r);
in.olabels    =    {'|a|^2','\int |da/dx|^2 dx'};
in.name =          'NLS soliton using finite differences + Dirichlet';
in.boundaries  =  [0,1;0,1];
in.da        =    @(a,w,r) 1i*a.*(conj(a).*a)+0.5*1i*(xd2(a,2,r)-a)+w*r.N;
e =xspde(in);                                           %%main program
end                                                     %%end of function