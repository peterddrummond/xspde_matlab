function [e] = SolitonDeriv()
%   e  =  SOLITONDERIV() tests xSPDE for a nonlinear Schroedinger equation.
%   Tests a (one+one)-dimensional partial differential equation for:
%   (1) using explicit derivatives in the derivative, in.da
%   (2) Using additional points and steps for higher accuracy
%   (3) Changing the integration method to RK4 with no linear term
%   (4) Adding an apodisation function to absorb at the boundaries
%   (5) Using derivatives to evaluate the knietic energy
%   (6) Using observables to plot an auxiliary function
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =          'NLS soliton using derivatives';
in.dimension =     2;                                     %%dimension: 1-4
in.points =        [101,81];                                     %%dimension: 1-4
in.steps =         40;
in.step =          @xRK4;
in.initial =       @(w,r)      2*sech(r.x);                 %%Initialisation
A        =         @(a,r)      (cos(pi*r.x/r.ranges(2)).^0.2-1); %%Apodisation 
in.da =            @(a,~,r)    i*a.*(conj(a).*a)+0.5*i*(xd(a,r.Dx.^2,r)-a)+a.*A(a,r);
in.observe{1} =    @(a,r)   a.*conj(a);
in.observe{2} =    @(a,r)   xint (abs(xd(a,r.Dx,r)).^2,r);
in.observe{3} =    @(a,r)   A(a,r);
in.olabels    =    {'a','\int [da/dx]^2 dx','A(x)'};
e =xspde(in);                                             %%main program
end                                                       %%end of function