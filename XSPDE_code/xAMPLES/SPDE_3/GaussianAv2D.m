function [e] = GaussianAv2D()
%   e  =  GAUSSIANAV2D() tests xSPDE for a linear Schroedinger equation
%   Tests a (two+one)-dimensional partial differential equation for:
%   (1) Initial gaussian in space
%   (2) Setting da to zero using the zeros function
%   (3) Averaging over the spatial grid
%   (4) Integrating over the spatial grid
%   (5) Integrating over the momentum grid
%   (6) Using transforms, images, transverse, compare, pdimension options
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name        =  'GaussianAv2D: diffraction in 2D';
p.dimensions  =  3;                                   %%dimension: 1-4
p.initial     =  @(~,p) exp(-0.5*(p.x.^2+p.y.^2));    %%initialisation
p.deriv       =  @(a,~,~) zeros(size(a));             %%da
p.observe{1}  =  @(a,p) a.*conj(a);                   %%observable 1
p.observe{2}  =  @(a,p) Ave(a.*conj(a),p.dx,p);       %%observable 2
p.observe{3}  =  @(a,p) Int(a.*conj(a),p.dx,p);       %%observable 3
p.observe{4}  =  @(a,p) a.*conj(a);                   %%observable 4
p.linear      =  @(p) 1i*0.05*(p.Dx.^2+p.Dy.^2);      %%laplacian
p.transforms  =  {0,0,0,[0,1,1]};                     %%transforms
p.images      =  {2};                                 %%number of images
p.transverse  =  {2};                                 %%transverse plots
p.olabels     =  {'I','<I>','\int I dV','I(k)'};      %%labels
p.compare{1}  =  @gaussian2d;                         %%comparison
p.compare{4}  =  @(r) 1;                              %%comparison
p.axes{4}     =  {0,-1,-1};                           %%comparison axes
p.transverse  =  {3,0,0,0};
p.pdimension  =  {3,1,1,1};
e             =  xspde(p);                            %%simulation
end                                                   %%end of main

function c = gaussian2d(p)
sig  =  (1+(p.t/10).^2);
r2   =  p.x.^2+p.y.^2;
c    =  sig.^(-1).*exp(-r2./sig);
end
