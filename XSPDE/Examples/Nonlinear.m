function e = Nonlinear()
%   e  =  NONLINEAR() tests xSPDE for a nonlinear SDE
%   Tests a one-dimensional nonlinear stochastic differential equation for:
%   (1) Using local and parallel, but no serial ensembles
%   (2) Integrating an SDE with the RK4 method
%   (3) An initial condition with a mean and a random term
%   (4) Calculating six different observables
%   (5) Plotting six observables
%   (6) Using xave to generate local averages
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'Nonlinear oscillator, RK4';
in.ensembles =  [2000,1,6];                  %%samples,ensembles,parallel
in.steps =      2;
in.points =     100;
in.step       = @xRK4;
in.initial   =  @(w,r)   0.5+0.25*w;        %%Initialisation  handle
in.da        =  @(a,xi,r) (a-a.^3)+xi;         %%Derivative  handle
in.observe    =  {@(a,r) xave(a),@(a,r) xave(a.^2), @(a,r) a.^3, ...
                 @(a,r) a.^4, @(a,r) a.^5, @(a,r) a.^6};
in.olabels    =  {'<a>','<a^2>','<a^3>','<a^4>','<a^5>','<a^6>'};
e            =  xspde(in);                  %%Stochastic program
end
