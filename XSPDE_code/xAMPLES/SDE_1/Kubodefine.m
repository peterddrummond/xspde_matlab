function [e] = Kubodefine()
%   e  =  KUBODEFINE() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Serial ensemble operation
%   (2) Setting additional steps to improve accuracy
%   (3) Using two exactly solvable observables
%   (4) Defining two comparison functions
%   (5) Changing the integration method to xMPdefine
%   (6) Changing the label for time to tau
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name        =  'Kubodefine: 3 observables, auxiliary fields';
p.auxfields   =  1;
p.ensembles   =  [100,10];                        %%samples,ensembles,parallel
p.initial     =  @(w,p)   1+0*w;                  %%Initialisation  handle
p.steps       =  5;
p.define      =  @(a,w,p) exp(5*1i*p.t).*a(1,:);
p.deriv       =  @(a,w,p) 1i*w.*a(1,:);           %%Derivative  handle
p.observe{1}  =  @(a,r) real(a(1,:));
p.transforms  =  {0,1,1};
p.observe{2}  =  @(a,~) a(1,:).*conj(a(1,:));
p.observe{3}  =  @(a,~) a(2,:).*conj(a(2,:));
p.olabels     =  {'a','|a(\omega)|^2','|FT(a*e^{(5it))}|^2'};
p.xlabels     =  {'\tau'};
p.compare{1}  =  @(p) exp(-p.t/2);                %%Comparison function
e             =  xspde(p);                        %%Stochasic program
end
