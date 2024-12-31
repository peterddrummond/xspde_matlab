function e = Nonlinear()
%   e  =  NONLINEAR() tests xSPDE for a nonlinear SDE
%   Tests a one-dimensional nonlinear stochastic differential equation for:
%   (1) Using local and parallel, but no serial ensembles
%   (2) Integrating an SDE with the RK4 method
%   (3) An initial condition with a mean and a random term
%   (4) Calculating six different observables
%   (5) Plotting six observables
%   (6) Using an inline cell array
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name       =  'Nonlinear oscillator, RK4';
p.ensembles  =  [2000,1,6];                  %%samples,ensembles,parallel
p.steps      =  2;
p.points     =  100;
p.method     =  @RK4;
p.initial    =  @(v,~)   0.5+0.25*v;         %%Initialisation  handle
p.deriv      =  @(a,w,r) (a-a.^3)+w;         %%Derivative  handle
p.observe    =  {@(a,p) a,   @(a,p) a.^2, @(a,p) a.^3, ...
		         @(a,p) a.^4, @(a,p) a.^5, @(a,p) a.^6};
p.olabels    =  {'<a>','< a^2 >','< a^3 >','< a^4 >','< a^5 >','< a^6 >'};
e            =  xspde(p);                    %%Stochastic program
end
