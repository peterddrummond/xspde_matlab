function [e] = GainTransfer()                      %%name of main function

in.name =      'Gain: Tests transfer function';%name for simulation
in.ranges =    4;                          %%ranges: t,x,y,z  
in.noises =    [2,0];                      %%xnoises, knoises per point
in.ensembles = [10000,4,1];                %%samples,ensembles,parallel
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation 
in.da =        @(a,xi,r)  -a + xi(1,:)+1i*xi(2,:);   %%Derivative  
in.observe =   {@(a,~) a.*conj(a)};        %%Observe  handle
in.olabels =   {'<|a|^2>'};                %%labels for observables
in.compare =   {@(t,~) 1+0*t};             %%Comparison handle
in2        =   in;                         %%Second input
in2.steps  =    2;                         %%Steps per plotted point
in2.name =     'Gain with noise';          %%name for simulation
in2.transfer =  @(w,r,a0,r0) a0/sqrt(2)+(w(1,:)+1i*w(2,:))/2;
in2.da =        @(a,xi,~)  a + xi(1,:)+1i*xi(2,:);   %%Derivative
in2.observe{2} =@(a,~) a;                  %%Observe  handle
in2.olabels{2} ='<a>';                     %%labels for observables
in2.compare =   {@(t,~) 2*exp(2*t)-1,@(t,~) 0*t}; %%Comparison handle
e           =   xspde({in,in2});           %%Stochastic program
end                                             
