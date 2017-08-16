function [e] = GainTransfer() 
%   e  =  GAINTRANSFER() tests xSPDE for a sequence of SDEs with transfer
%   Tests a one-dimensional linear stochastic differential equation for:
%   (1) Transferring an attenuated field with additional noise
%   (2) Reusing the first observe function and label
%   (3) Introducing a second observe function for the second simulation
%   (4) Setting one compare for the first sequence, two for the second
%   (5) Testing the RK2 method for the first sequence
%   (6) Testing the RK4 method for the second sequence
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =      'Gain: Tests transfer function';%name for simulation
in.ranges =    4;                          %%ranges: t,x,y,z 
in.step   =    @xRK2;                      %%first integration method
in.noises =    [2,0];                      %%xnoises, knoises per point
in.ensembles = [1000,4,1];                 %%samples,ensembles,parallel
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2); %%Initialisation 
in.da =        @(a,xi,r)  -a + xi(1,:)+1i*xi(2,:);   %%Derivative  
in.observe =   {@(a,~) a.*conj(a)};        %%Observe  handle
in.olabels =   {'<|a|^2>'};                %%labels for observables
in.compare =   {@(t,~) 1+0*t};             %%Comparison handle
in2        =   in;                         %%Second input
in2.step   =   @xRK4;                      %%Second integration method
in2.steps  =    2;                         %%Steps per plotted point
in2.name =     'Gain with noise';          %%name for simulation
in2.transfer =  @transnoise;
in2.da =        @(a,xi,~)  a + xi(1,:)+1i*xi(2,:);   %%Derivative
in2.observe{2} =@(a,~) real(a);                  %%Observe  handle
in2.olabels{2} ='<a>';                     %%labels for observables
in2.compare =   {@(t,~) 2*exp(2*t)-1,@(t,~) 0*t}; %%Comparison handle
e           =   xspde({in,in2});           %%Stochastic program
end  

function [a0,r] = transnoise(w,r,a,~)       %% Subsequent data
    a0 = a/sqrt(2)+(w(1,:)+1i*w(2,:))/2;    %% set to previous output
end
