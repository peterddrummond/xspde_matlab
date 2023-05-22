function [e] = GaussianAv4D()
%   [e]  =  GAUSSIANAV4D() tests xSPDE for a linear Schroedinger equation
%   Tests a (four+one)-dimensional partial differential equation for:
%   (1) Initial 4D gaussian in space, using numeric grid axis notation
%   (2) Setting da to zero using the default
%   (3) Averaging over the spatial grid using one argument for xave
%   (4) Integrating over the spatial grid using two arguments for xint
%   (5) Averaging and integrating over the momentum grid
%   (6) Using transforms, images, transverse, compare, pdimension in 4+1 D
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'GaussianAv4D: 4+1D Schroedinger equation';
p.dimensions   =  5;
p.noises       =  2;
p.points       =  [7,15,15,15,15];
p.initial      =  @(w,p) exp(-0.5*(p.r{2}.^2+p.r{3}.^2+p.r{4}.^2+p.r{5}.^2));
p.observe{1}   =  @(a,p) a.*conj(a);                            %%observable 1
p.observe{2}   =  @(a,p) Ave(a.*conj(a),p);                     %%observable 2
p.observe{3}   =  @(a,p) Int(a.*conj(a),p);                     %%observable 3
p.observe{4}   =  @(a,p) a.*conj(a);                            %%observable 4
p.observe{5}   =  @(a,p) Ave(a.*conj(a),p);                     %%observable 5
p.observe{6}   =  @(a,p) Int(a.*conj(a),p.dk,p);                %%observable 6
p.linear       =  @(p) 1i*0.05*(p.D{2}.^2+p.D{3}.^2+p.D{4}.^2+p.D{5}.^2);
p.transforms   =  {0,0,0,[0,1,1,1,1],[0,1,1,1,1],[0,1,1,1,1]};
p.images       =  2;                                            %%images
p.transverse   =  {2,0,0,2};                                    %%transverse plots
p.olabels      =  {'I','<I>','\int I dV','I(k)','<I(k)>','\int I dK'};%%label
p.compare{1}   =  @gaussian4d;                                  %%comparison
p.compare{4}   =  @(p) 1;                                       %%comparison
p.axes{4}      =  {0,-1,-1,-1,-1};                              %%comp. axes
p.pdimension   =  {3,1,1,2,1,1};
e              =  xspde(p);
end                                                             %%end of main


function c = gaussian4d(p)
sig  =  (1+(p.t/10).^2);
r2   =  p.r{2}.^2+p.r{3}.^2+p.r{4}.^2+p.r{5}.^2;
c    =  sig.^(-2).*exp(-r2./sig);
end
