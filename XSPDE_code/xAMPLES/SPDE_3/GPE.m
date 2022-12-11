function [e] = GPE()
%   e  =  GPE() tests xSPDE for  a vortex in two space dimensions
%   Tests a three-dimensional noisy GPE in a harmonic well by
%   (1) Modifying points and ranges in 3 dimensions
%   (2) Using momentum filters on noises generated
%   (3) Setting the noise and random transform filter in k-space
%   (4) Changing the movie image type to give a contour plot
%   (5) Using auxiliary functions defined inline
%   (6) Graphing data from results saved in an HDF5 file
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name        =  'GPE';
p.dimensions  =  3;
p.points      =  [101,63,63];
p.ranges      =  [1,20,20];
p.noises      =  [0,2];
p.rfilter     =  @(w,p)   w.*exp(- p.kx.^2-p.ky.^2);
p.nfilter     =  @(v,p)   v.*exp(- p.kx.^2-p.ky.^2);
b             =  @(w,p)   reshape(.1*(w(1,:)+1i*w(2,:)),p.d.r);
p.initial     =  @(w,p)   (p.x+1i*p.y)./(1+10*(p.x.^2 +p.y.^2))+b(w,p);
V             =  @(p)     0.01*(p.x.^2 +p.y.^2)-0.001*1i*(p.x.^2 +p.y.^2).^2;
p.deriv       =  @(a,w,p) -1i*a.*(V(p)+conj(a).*a)+b(w,p);
p.linear      =  @(p)     0.5*1i*(p.Dx.^2+p.Dy.^2);
p.observe{1}  =  @(a,p)   a.*conj(a);
p.images      =  {2};
p.imagetype   =  {2};
p.olabels     =  {'|a|^2'};
p.file        =  'GPE.mat';
[e,~,p]       =  xsim(p);                  %%sim program
xgraph(p.file,p);                          %%graph program
end                                        %%end of function


