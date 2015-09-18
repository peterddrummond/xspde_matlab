function [e] = Kubotest()                         

in.name =       'Kubotest: 2 observables';
in.ensembles =  [2000,10];                       %%samples,ensembles,parallel
in.initial   =  @(w,r)    1+0*w ;              %%Initialisation  handle
in.steps     =  5;
in.da        =  @(a,xi,r) 1i*xi.*a;             %%Derivative  handle
in.observe{1} = @(a,r) a(1,:);
in.observe{2} = @(a,r) a(1,:).*conj(a(1,:));
in.olabels    = {'a','|a|^2'};
in.compare{1} = @(t,~) exp(-t/2);               %%Comparison function
in.compare{2} = @(t,~) 1+0*t;                       %%Comparison function
[e,data]  =     xspde(in);                   %%Stochasic program
end
