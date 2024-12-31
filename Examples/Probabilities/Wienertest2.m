function e = Wienertest2()
%   e  =  WIENERTEST2() tests xSPDE for a basic Wiener process
%   Tests a two-variable additive SDE equation for:
%   (1) Binning for a 2D probability with a different scatter plot
%   (2) Transverse probability plots with two variables
%   (3) Error comparison with a compare function in two dimensions
%   (4) Using a displaced variable for the probability argument
%   (5) Multivariate labels with subscripts
%   (6) Multiple 2D and 3D probability plots with different labels and axes
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name           =  'Wienertest2';                  %%name of simulation
p.method         =  @Euler;                         %%Euler method
p.fields         =  2;                              %%two components
p.noises         =  2;                              %%noises set to two
p.points         =  21;                             %%time points
p.ensembles      =  [10000,10];                      %%ensembles 
p.initial        =  @(v,p) 0.5*v;                   %%Initial function
p.deriv          =  @(a,w,p) [0.5*w(1,:);2*w(2,:)]; %%Derivative function
p.observe{1}     =  @(a,~) a;                       %%Observable function
p.binranges{1}   =  {-5:0.25:5,-5:0.25:5};          %%Bin boundaries
p.cutoffs{1}     =  .004;                           %%Chi-square cutoffs
p.transverse{1}  =  3;                              %%Transverse plots
p.images{1}      =  3;                              %%Image plots
p.glabels{1}     =  {'t','a_1','a_2'};              %%Probability labels
p.compare{1}     =  @wienercp2;                     %%Observable function
p.scatters{2}    =  12;                             %%Scatters function
p.observe{2}     =  @(a,~) a;                       %%Observable function
p.olabels        =  {'P(a)','a_1'};                 %%vertical axis labels
e                =  xspde(p);                       %%Runs xspde
end                                                 %%end of main function

function c = wienercp2(p)
sig1  =  0.25*(1+p.t);
sig2  =  0.25+4.*p.t;
nm    =  (4*pi^2*sig1.*sig2).^(-0.5);
c     =  nm.*exp(-(p.r{2}.^2./(2*sig1)));
c     =  c.*exp(-(p.r{3}.^2./(2*sig2)));
end
