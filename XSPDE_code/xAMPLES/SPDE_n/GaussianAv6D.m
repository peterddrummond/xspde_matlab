function [e]  =  GaussianAv6D()
%   [e]  =  GAUSSIANAV6D() tests xSPDE for a linear Schroedinger equation
%   Tests a (six+one)-dimensional partial differential equation for:
%   (1) Initial 6D gaussian in space, using numeric grid axis notation
%   (2) Setting da to zero using the default in 6D
%   (3) Integrating over part of the spatial grid using xint
%   (5) Integrating over part of the the momentum grid
%   (6) Using transforms over part of the spatial grid
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name          =  '6+1D Schroedinger equation';
p.dimensions    =  7;
p.points        =  [2,9,9,9,9,9,9];                           %%dimension: 1-4
p.ranges        =  [10,7,7,7,7,7,7];
p.initial       =  @(v,p) exp(-0.5*(p.r{2}.^2+p.r{3}.^2+p.r{4}.^2+p.r{5}.^2+p.r{6}.^2+p.r{7}.^2));
p.observe{1}    =  @(a,p) a.*conj(a);
p.observe{2}    =  @(a,p) Int(a.*conj(a),p.dx.*[0,1,1,0,0,1,1],p);
p.observe{3}    =  @(a,p) a.*conj(a);
p.observe{4}    =  @(a,p) Int(a.*conj(a),[p.dx(1:3),p.dk(4:7)],p);
p.linear        =  @(p) 1i*0.05*(p.D{2}.^2+p.D{3}.^2+p.D{4}.^2+p.D{5}.^2+p.D{6}.^2+p.D{7}.^2);
p.transforms    =  {0,0,[0,0,0,1,1,1,1],[0,0,0,1,1,1,1]};
p.images        =  3;                                         %%number of images
p.imagetype     =  2;                                         %%number of images
p.transverse    =  1;                                         %%transverse plots
p.raw           =  1;
%in.file        =  'GaussianAv6D.mat';   %%Warning - makes large disk file!
p.olabels       =  {'I','\int I dx_1dx_2dx_5dx_6','I(k)','\int I dK'};
p.axes{1}       =  {-1,-1,-1,-1,0,0,0};
p.axes{2}       =  {0,-1,-1,-1,-1,-1,-1};
p.axes{3}       =  {0,0,0,-1,-1,-1,-1};
p.compare{2}    =  @(p) pi^2*(1+(p.t/10).^2).^(-1);           %%comparison
%in.compare{3}  =  @(x) (1+(x.t/10).^2).^(-1);                %%comparison
p.compare{3}    =  @gaussian6dft;
p.compare{4}    =  @(p) pi^3;                                 %%comparison
p.pdimension    =  {3,1,1,1,1,1};
e               =  xspde(p);                                  %%simulation
end                                                           %%end of main

function c = gaussian6dft(p)
sig  =  (1+(p.t/10).^2);
r2   =  p.r{2}.^2+p.r{3}.^2;
c    =  sig.^(-1).*exp(-r2./sig);
end
