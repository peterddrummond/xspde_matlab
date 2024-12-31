function [e] = Kubo2()
%   e  =  KUBO2() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Serial ensemble operation
%   (2) Setting additional steps to improve accuracy
%   (3) Using two exactly solvable observables
%   (4) Defining two comparison functions
%   (5) Changing the integration method to xRK2
%   (6) Changing the label for time to tau
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'Kubo2: 2 observables';
p.ensembles    =  [1000,10];                       %%samples,ensembles,parallel
p.initial      =  @(w,p)    1+0*w;                 %%Initialisation  handle
p.method       =  @RK2;
p.steps        =  5;
p.order        =  1;
p.deriv        =  @(a,w,~) 1i*w.*a;                %%Derivative  handle
p.observe{1}   =  @(a,~) real(a(1,:));
p.observe{2}   =  @(a,~) a(1,:).*conj(a(1,:));
p.olabels      =  {'a','|a|^2'};
p.xlabels      =  {'\tau'};
p.compare{1}   =  @(p) exp(-p.t/2);                %%Comparison function
p.compare{2}   =  @(p) 1+0*p.t;                    %%Comparison function
e              =  xspde(p);                        %%Stochasic program
end
