function [e] = Kubotest()                            
%   e  =  KUBO() tests xSPDE for a Kubo oscillator
%   Tests a one-variable multiplicative SDE equation for:
%   (1) Serial ensemble operation
%   (2) Setting additional steps to improve accuracy
%   (3) Using two exactly solvable observables
%   (4) Defining two comparison functions
%   (5) Changing the integration method to xRK2
%   (6) Changing the label for time to tau
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License                         

in.name =       'Kubotest: 2 observables';
in.ensembles =  [2000,10];                       %%samples,ensembles,parallel
in.initial   =  @(w,r)    1+0*w ;                %%Initialisation  handle
in.step      =  @xRK2;
in.steps     =  5;
in.da        =  @(a,xi,r) 1i*xi.*a;              %%Derivative  handle
in.observe{1} = @(a,r) a(1,:);
in.observe{2} = @(a,r) a(1,:).*conj(a(1,:));
in.olabels    = {'a','|a|^2'};
in.xlabels =    {'\tau'};
in.compare{1} = @(t,~) exp(-t/2);                %%Comparison function
in.compare{2} = @(t,~) 1+0*t;                    %%Comparison function
e         =     xspde(in);                       %%Stochasic program
end
