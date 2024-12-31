function [e] = Kubodefine()
%   e  =  KUBODEFINE() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Serial ensemble operation
%   (2) Setting additional steps to improve accuracy
%   (3) Using an auxiliary field observable function
%   (4) Defining a comparison functions
%   (5) Using p.verbose = 1 for detailed outputs
%   (6) Changing the label for time to tau
%   Licensed by Peter D. Drummond, (2024) - see License

p.name        =  'Kubodefine';
p.auxfields   =  1;
p.ensembles   =  [100,10];                       %%samples,ensembles,par.
p.initial     =  @(w,p)   1+0*w;                 %%Initialisation  handle
p.steps       =  5;
p.verbose     =  1;
p.checks = 1
p.define      =  @(a,w,p) exp(5*1i*p.t).*a;
p.deriv       =  @(a,w,p) 1i*w.*a;               %%Derivative  handle
p.observe{1}  =  @(a,x,r) a;
p.transforms  =  {0,1,1};
p.observe{2}  =  @(a,x,~) a.*conj(a);
p.observe{3}  =  @(a,x,~) x.*conj(x);
p.olabels     =  {'a','|a(\omega)|^2','|FT(a*e^{(5it))}|^2'};
p.xlabels     =  {'\tau'};
p.compare{1}  =  @(p) exp(-p.t/2);                %%Comparison function
e             =  xspde(p);                        %%Stochasic program
end
