function e         =  Wienertest()
%   e  =  WIENERTEST() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for:
%   (1) Binning to get a probability
%   (2) Transverse probability plots
%   (3) Error comparison with a compare function
%   (4) Observed variable for the probability argument
%   (5) Binning to get a second probability observable
%   (6) Multiple 2D and 3D probability plots with different labels
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name            =  'Wienertest';             %%name of simulation
p.method          =  @Euler;                   %%first integration method
p.noises          =  1;                        %%noises set to one
p.ensembles       =  [1000,10];                %%ensembles for averaging
p.initial         =  @(v,p)   v;               %%Initial function
p.deriv           =  @(a,w,p) w;               %%Derivative function
p.observe{1}      =  @(a,~) a;                 %%Observable function
p.binranges{1}    =  {-5:0.25:5};              %%Bin boundaries
p.transverse{1}   =  2;                        %%Transverse plots
p.cutoffs         =  {0.004,0.002};            %%Cutoffs
p.compare{1}      =  @wienerc1;                %%Compare
p.olabels         =  {'P(a)','P(a^2)'};        %%labels
p.glabels         =  {{'t','a'},{'t','a^2'}};  %%labels
p.observe{2}      =  @(a,~) a.^2;              %%Observable function
p.binranges{2}    =  {0.5:0.5:25};             %%Bin boundaries
p.transverse{2}   =  2;                        %%Transverse plots
p.compare{2}      =  @wienerc2;                %%Compare
e                 =  xspde(p);                 %%Runs xspde simulation
end                                            %%end of main function

function c1 = wienerc1(p)
sig  =  1+p.t;
nm   =  (2*pi*sig).^(-0.5);
c1   =  nm.*exp(-(p.r{2}.^2./(2*sig)));
end

function c2 = wienerc2(p)
sig  =  1+p.t;
nm   =  (2*pi*sig).^(-0.5);
c2   =  (nm./sqrt(p.r{2})).*exp(-(p.r{2}./(2*sig)));
end
