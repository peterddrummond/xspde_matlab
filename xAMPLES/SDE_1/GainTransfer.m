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

p.name          =  'Gain Transfer I';                   %%name for simulation
p.ranges        =  4;                                   %%ranges: t,x,y,z
p.method        =  @RK2;                                %%first integration method
p.noises        =  2;                                   %%xnoises, knoises per point
p.ensembles     =  [500,1,12];                          %%samples,ensembles,parallel
p.initial       =  @(w,~) (w(1,:)+1i*w(2,:))/sqrt(2);   %%Initialisation
p.deriv         =  @(a,w,p)  -a + w(1,:)+1i*w(2,:);     %%Derivative
p.observe       =  {@(a,~) a.*conj(a)};                 %%Observe  handle
p.olabels       =  {'<|a|^2>'};                         %%labels for observables
p.compare       =  {@(r) 1};                            %%Comparison handle
p2              =  p;                                   %%Second input
p2.method       =  @RK4;                                %%Second integration method
p2.steps        =  2;                                   %%Steps per plotted point
p2.name         =  'Gain Transfer II';                   %%name for simulation
p2.transfer     =  @transnoise;
p2.deriv        =  @(a,w,~)  a + w(1,:)+1i*w(2,:);      %%Derivative
p2.observe{2}   =  @(a,~) real(a);                      %%Observe  handle
p2.olabels{2}   =  '<a>';                               %%labels for observables
p2.compare      =  {@(p) 2*exp(2*(p.t-p.origins(1)))-1,@(p) 0}; %%Comparison
e               =  xspde(p,p2);                       %%Stochastic program
end

function [a0,p] = transnoise(a,w,p)                     %% Subsequent data
    a0 = a/sqrt(2)+(w(1,:)+1i*w(2,:))/2;                %% set to previous output
end
