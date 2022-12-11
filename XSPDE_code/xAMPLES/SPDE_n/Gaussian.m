function [e]  =  Gaussian()
%   [e]  =  GAUSSIAN() tests xSPDE for a linear Schroedinger equation
%   Tests a four-dimensional partial differential equation for:
%   (1) Inputting the x,y,z-grid into the initial condition
%   (2) Default derivative function set to zero in a PDE
%   (3) Setting up the 3D linear response including a linear loss term
%   (4) Inputting two compare functions and calculating differences
%   (5) Graphing results with transverse plots
%   (6) Graphing results with movie images
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name         =  'Gaussian diffraction with loss';
p.dimensions   =  4;                                            %%dimension: 1-4
p.initial      =  @(w,p) exp(-(p.x.^2+p.y.^2+p.z.^2)/2);        %%initialisation
p.observe{1}   =  @(a,~) a.*conj(a);                            %%x-observable
p.linear       =  @(p) -0.5+1i*0.05*(p.Dx.^2+p.Dy.^2+p.Dz.^2);  %%laplacian
p.images       =  {2,2};                                        %%number of images
p.transverse   =  {2,2};                                        %%transverse plots
p.olabels      =  {'|a(r)|^2'};                                 %%labels
p.compare{1}   =  @gaussian3d;                                  %%comparison
e              =  xspde(p);                                     %%simulation
end                                                             %%end of main

function c  =  gaussian3d(p)
sig         =  (1+(p.t/10).^2);
r2          =  p.x.^2+p.y.^2+p.z.^2;
c           =  sig.^(-3/2).*exp(-p.t-r2./sig);
end
