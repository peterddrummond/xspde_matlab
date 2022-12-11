function [e]        =  SHO()
%   [e]  =  SHO() tests xSPDE for a simple harmonic oscillator
%   Tests a one-variable differential equation for:
%   (1) Vector calculations
%   (2) Exactly soluble DE case
%   (3) Default range, points, dimension
%   (4) Adding an initial value
%   (5) Default observe functions, olabels
%   (6) Comparisons, scatters,parametric plots, sequential plots
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 p.name             =  'Simple harmonic oscillator';        %%Name
 p.initial          =  @(w,~) [1;0];                        %%Initialisation
 p.fields           =  2;
 p.steps            =  20;
 p.order            =  2;
 p.ensembles        =  3;
 p.headers          =  {'x-plot','y-plot','parametric plot'};
 p.olabels          =  {'x','y','x'};
 p.deriv            =  @(a,w,~) [-a(2,:);a(1,:)];           %%Derivative function
 p.compare          =  {@(p) cos(p.t)};                     %%Compare  handle
 p.observe{1}       =  @(a,p) a(1,:);
 p2                 =  p;
 p2.deriv           =  @(a,w,~) [-a(2,:);a(1,:)]+0.005*w;   %%Derivative function
 p2.compare{1}      =  [];
 p2.observe{2}      =  @(a,~) a(2,:);
 p2.observe{3}      =  @(a,~) a(1,:);
 p2.limit{1}        =  [5,10];
 p2.scatters        =  {1,3,3};
 p2.parametric{3}   =  [2,1];
 e                  =  xspde({p,p2});
end                                                         %%end of main function
