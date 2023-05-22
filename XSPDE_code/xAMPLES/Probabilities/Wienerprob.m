function e = Wienerprob()
%   e  =  WIENERPROB() tests xSPDE for a basic Wiener process
%   Tests a one-variable additive SDE equation for:
%   (1) Plotting multiple lines on a graph
%   (2) Using legends on an output graph
%   (3) Setting up the probability bin ranges
%   (4) Multiple line comparisons with exact results
%   (5) Using probabilities with comparisons and chi-squares
%   (6) Defining an external comparison probability function
%   xSPDE is licensed by Peter D. Drummond, (2022) - see License

p.name           =  'Wienerprob: SDE distribution';  %%name of simulation
p.noises         =  1;                               %%noises set to one
p.points         =  10;                              %%number of time points
p.ensembles      =  [10000,10];                      %%ensembles for averaging
p.initial        =  @(v,p) v/2;                      %%Initial function
p.sig            =  @(p) .25 + p.r{1};               %%variance
p.deriv          =  @(a,w,p) w;                      %%Derivative function
p.observe{1}     =  @(a,p) a;                        %%Observable function
p.compare{1}     =  @gaussprob;                      %%Comparison function
p.observe{2}     =  @(a,p) [a;a.^2;a.^3;a.^4/30];    %%Observable function
p.compare{2}     =  @(p) [0*p.t;.25+p.t;0*p.t;.1*(.25+p.t).^2]; %%Exact
p.transverse{1}  =  5;                               %%Transverse graphs
p.olabels{1}     =  'P(x)';                          %%function labels
p.binranges{1}   =  {-5:0.25:5};                     %%binning ranges
p.legends{1}     =  {'Sampled P(x,\tau) \pm \sigma','Exact  P(x,\tau)'};
p.olabels{2}     =  '<x^n>';
p.xlabels        =  {'\tau','x'};
p.legends{2}     =  {'<x>','<x^2>','<x^3>','<x^4>'};
e                =  xspde(p);                        %%Runs xspde simulation
end                                                  %%end of main function

function p = gaussprob(p)
p  =  exp(-(p.r{2}.^2)./(2*p.sig(p)))./sqrt(2*pi*p.sig(p));
end
